nextflow.enable.dsl=2

//to run must have nextflow and miniconda3 loaded.
//see yaml of conda env packages included in project
//also need to install AlphaPlantImpute2 separately 

//  Set up parameters
// tools
params.conda = '/home/nicolas.lara/.conda/envs/imputation'
// directories
params.vcf_dir = '/project/guedira_seq_map/nico/vcf_files'
params.output_dir = '/90daydata/guedira_seq_map/plink_output'
// keyfiles
params.chrom_remap = '/project/guedira_seq_map/nico/chromosome_remap.txt'
params.parents = '/project/guedira_seq_map/nico/sunrils_parents.txt'


process process_vcfs {
//    conda '${params.conda}'
    conda '/home/nicolas.lara/.conda/envs/imputation'

    input:
    path vcf

    output:
    tuple path("${sample}_parents.vcf"), path("${sample}_pop.vcf")

    script:
    sample = vcf.baseName
    """
    ##rename chromosomes to numbers for downstream processing
    bcftools annotate --rename-chrs ${params.chrom_remap} ${vcf} > vcf1.vcf
    ##remove unknown/unaligned markers, filter to biallelic with MAF > 0.05
    bcftools view -t "^UN,Unknown" vcf1.vcf > vcf2.vcf
    bcftools view -m2 -M2 -v snps -i'MAF>0.05' vcf2.vcf > vcf3.vcf
    ##filter to one chromosome for testing
    bcftools view -t 1 vcf3.vcf > vcf4.vcf
    mv vcf4.vcf vcf3.vcf
    ##split into parent and population files
    bcftools view -S ${params.parents} vcf3.vcf > ${sample}_parents.vcf
    bcftools view -S ^${params.parents} vcf3.vcf > ${sample}_pop.vcf
    ##clean up
    rm vcf*.vcf
    """
}

process mask_pop_vcfs {
    input:
    path vcf_file

    output:
    tuple path("pop_masked.vcf.gz"), path("pop_antimasked.vcf.gz")

    conda '${CONDA_ENV}'

    script:
    sample = vcf.baseName
    """
    #!/usr/bin/env Rscript
    library(vcfR)
    ##read in file and extract markers
    vcf <- read.vcfR("${vcf_file}")
    gt <- extract.gt(vcf, as.numeric=T)
    ##split into two matrices, masking a complementary n% of each
    n <- 0.5
    masks <- apply(gt, MARGIN = 2, function(x) {
        na_indices <- sample(length(x), length(x) * n)
        x_inverse <- x
        x[na_indices] <- NA
        x_inverse[-na_indices] <- NA
        list(mask = x, antimask = x_inverse)
    }
    ##turn into two matrices
    mask_m <- do.call(cbind, lapply(masks, function(col), col$mask)
    anti_m <- do.call(cbind, lapply(masks, function(col), col$antimask)
    colnames(mask_m) <- colnames(gt); rownames(mask_m) <- rownames(gt)
    colnames(anti_m) <- colnames(gt); rownames(anti_m) <- rownames(gt)
    ##coerce back into vcfR object and save
    mask_vcf <- vcf; anti_vcf <- vcf
    mask_vcf@fix <- mask_vcf@fix[mask_vcf@fix[,"ID"] %in% rownames(mask_m),]
    mask_vcf@gt <- cbind(as.matrix(data.frame(FORMAT = rep("GT", nrow(mask_m)))), mask_m)
    write.vcf(mask_vcf, file="pop_masked.vcf.gz")
    anti_vcf@fix <- anti_vcf@fix[anti_vcf@fix[,"ID"] %in% rownames(anti_m),]
    anti_vcf@gt <- cbind(as.matrix(data.frame(FORMAT = rep("GT", nrow(anti_m)))), anti_m)
    write.vcf(anti_vcf, file="pop_antimasked.vcf.gz")
    """
}

process vcf_to_plink {
    conda '/home/nicolas.lara/.conda/envs/imputation'

    input:
    tuple path(parent_vcf), path(population_vcf)

    output:
    tuple path("${parent_vcf.baseName}.ped"), path("${population_vcf.baseName}.ped")


    script:
    """
    plink --vcf ${parent_vcf} --recode --keep-allele-order --out ${parent_vcf.baseName}
    plink --vcf ${population_vcf} --recode --keep-allele-order --out ${population_vcf.baseName}
    """
}

process impute {
    publishDir "${params.output_dir}", mode: 'copy'
    conda '/home/nicolas.lara/.conda/envs/imputation'

    input:
    tuple path(parent_plink), path(population_plink)

    output:
    path("${population_plink.baseName}_imputed*"), emit: imputed


    script:
    """
    AlphaPlantImpute2 -createlib -out library -ped ${parent_plink}
    AlphaPlantImpute2 -impute -out ${population_plink.baseName}_imputed -library library.ped -ped ${population_plink}
    """
}

process plink_to_vcf {
    publishDir "{params.output_dir}", mode: 'copy'
    input:
    path plink

    output:
    path vcf

    conda '${CONDA_ENV}'

    script:
    sample = plink.baseName
    """
    plink --bfile ${plink} --recode vcf --out ${sample}
    """
}

process accuracy_evaluation {
    input:
    path masked_vcf
    path truth_vcf

    output:
    path "*.txt", emit accuracy

    conda '${CONDA_ENV}'

    script:
    """
    vcf-compare ${masked_vcf} ${truth_vcf} > results.txt
    bcftools isec ${truth_vcf} ${masked_vcf} -n =2 -p isec_output/
    """

}

workflow {
    Channel.fromPath(params.vcf_dir + '/*.vcf') \
        .set { vcf_files }

   //also add code to mask vcf and pass through.
   //so structure should be sample1: parents, pop, masked_pop
    processed_vcfs = process_vcfs(vcf_files) |
	vcf_to_plink

    plink_result = processed_vcfs |
	impute

    //take the triple files, create library for each sample from parents
    //then impute both the pop and masked_pop
    //imputed = impute(processed_vcfs)

    //take the masked_pops and turn back into vcfs
    //imputed_vcf = plink_to_vcf(imputed)
    
    //match imputed masked vcfs to raw population vcfs, compare
    //accuracy_evaluation(imputed_vcf)
}

