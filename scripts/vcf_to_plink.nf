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
    bcftools annotate --rename-chrs ${params.chrom_remap} ${vcf} > vcf1.vcf
    bcftools view -t "^UN,Unknown" vcf1.vcf > vcf2.vcf
    bcftools view -m2 -M2 -v snps -i'MAF>0.05' vcf2.vcf > vcf3.vcf
    bcftools view -t 1 vcf3.vcf > vcf4.vcf
    mv vcf4.vcf vcf3.vcf
    bcftools view -S ${params.parents} vcf3.vcf > ${sample}_parents.vcf
    bcftools view -S ^${params.parents} vcf3.vcf > ${sample}_pop.vcf
    rm vcf*.vcf
    """
}

process mask_pop_vcfs {
    input:
    path vcf

    output:
    path '*.vcf', emit: vcf

    conda '${CONDA_ENV}'

    script:
    sample = vcf.baseName
    """
    import random
    import pysam

    # Open the VCF file
    vcf = pysam.VariantFile("truth_data.vcf", "r")  # Open the VCF in read mode
    output_vcf = pysam.VariantFile("masked_90.vcf", "w", header=vcf.header)

    # Loop through records and randomly mask genotypes
    for rec in vcf.fetch():
        for sample in rec.samples:
            if random.random() < 0.9:  # Mask 90% of the genotypes
                rec.samples[sample]["GT"] = (None,)  # Set genotype to missing
        output_vcf.write(rec)

    output_vcf.close()
    vcf.close()
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
    plink --vcf ${parent_vcf} --recode --out ${parent_vcf.baseName}
    plink --vcf ${population_vcf} --recode --out ${population_vcf.baseName}
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

