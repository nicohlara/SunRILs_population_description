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
    tuple path(vcf_parents), path(vcf_pop)

    output:
    tuple path(vcf_parents), path("${sample}_masked.vcf.gz"), path("${sample}_antimasked.vcf.gz")

    conda '${CONDA_ENV}'

    script:
    sample = vcf_pop.baseName
    """
    #!/usr/bin/env Rscript
    library(vcfR)
    ##read in file and extract markers
    vcf <- read.vcfR("${vcf_pop}")
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
    write.vcf(mask_vcf, file="${sample}_masked.vcf.gz")
    anti_vcf@fix <- anti_vcf@fix[anti_vcf@fix[,"ID"] %in% rownames(anti_m),]
    anti_vcf@gt <- cbind(as.matrix(data.frame(FORMAT = rep("GT", nrow(anti_m)))), anti_m)
    write.vcf(anti_vcf, file="${sample}_antimasked.vcf.gz")
    """
}

process vcf_to_plink {
    conda '/home/nicolas.lara/.conda/envs/imputation'

    input:
    tuple path(parent_vcf), path(pop_vcf_mask), path(pop_vcf_anti)

    output:
    tuple path("${parent_vcf.baseName}.ped"), path("${pop_vcf_mask.baseName}.ped"), path(pop_vcf_anti)


    script:
    """
    plink --vcf ${parent_vcf} --recode --keep-allele-order --out ${parent_vcf.baseName}
    plink --vcf ${pop_vcf_mask} --recode --keep-allele-order --out ${pop_vcf_mask.baseName}
    """
}

process impute {
    publishDir "${params.output_dir}", mode: 'copy'
    conda '/home/nicolas.lara/.conda/envs/imputation'

    input:
    tuple path(parent_plink), path(population_plink), path(pop_vcf_anti)

    output:
    tuple path("${population_plink.baseName}_imputed*"), path(pop_vcf_anti)


    script:
    """
    AlphaPlantImpute2 -createlib -out library -ped ${parent_plink}
    AlphaPlantImpute2 -impute -out ${population_plink.baseName}_imputed -library library.ped -ped ${population_plink}
    """
}

process plink_to_vcf {
    publishDir "{params.output_dir}", mode: 'copy'
    input:
    tuple path(plink), path(pop_vcf_anti)

    output:
    tuple path(${sample}.vcf), path(pop_vcf_anti)

    conda '${CONDA_ENV}'

    script:
    sample = plink.baseName
    """
    plink --bfile ${plink} --recode vcf --out ${sample}
    """
}

process check_imp_accuracy {
    input:
    tuple path(imp_vcf), path(anti_vcf)

    output:
    path ${sample}_impute_accuracy.txt

    conda '${CONDA_ENV}'

    script:
    sample = vcf_pop.baseName
    """
    #!/usr/bin/env Rscript
    library(vcfR)
    ##read in file and extract markers
    imputed <- read.vcfR(${imp_vcf})
    complement <- read.vcfR(${anti_vcf})
    ##for GBS data
    gti <- extract.gt(imputed)
    gtc <- extract.gt(complement, as.numeric = T)    
    markers <- complement@fix[, c("ID", "REF", "ALT")]
    markers_good <- markers[markers[, "REF"] == imputed@fix[, "REF"] &
    markers[, "ALT"] == imputed@fix[, "ALT"], ]
    markers_flipped <- markers[markers[, "REF"] == imputed@fix[, "ALT"] &
    markers[, "ALT"] == imputed@fix[, "REF"], ]

    ## recode into major/hetero/minor
    gti[gti == "1/1"] <- 2
    gti[gti == "0/1"] <- 1
    gti[gti == "0/0"] <- 0
    mode(gti) <- "numeric"
    gtc[gtc == "1"] <- 2
    gtc[gtc == "0"] <- 0

    ## flip markers that need it
    gti[rownames(gti) %in% markers_flipped[, "ID"], ] <- abs(gti[rownames(gti) %in% markers_flipped[, "ID"], ] - 2)

    ## compare
    gti_align <- gti[rownames(gtc), colnames(gtc)]
    check_positions <- !is.na(gtc) & !is.na(gti_align)
    combination_table <- table(gtc[check_positions], gti_align[check_positions])
    combination_table
    total <- sum(combination_table)
    per_comb <- round(combination_table / total, 3)
    output <- as.vector(per_comb)
    names(output) <- c("major/major", "minor/major", "major/het", "minor/het", "major/minor", "minor/minor")
    write.table(data.frame(t(a)), file="${sample}_accuracy.csv", quote=F, sep=",", row.names=F)
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

