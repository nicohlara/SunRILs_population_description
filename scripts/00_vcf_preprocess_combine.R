nextflow.enable.dsl=2

params.vcf_files = "data/*.vcf.gz"
params.output_dir = "results/merged"
params.chromosomes = ["1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"]
// bcftools filtering parameters
params.depth = '3'
params.quality = '20'

workflow {
    // Create output directory
    new File(params.output_dir).mkdirs()
	
    Channel.fromPath('path/to/*.vcf.gz')
        .set { raw_vcfs }

    // Step 1: preprocess
    preprocessed = bcftools_filter(raw_vcfs)

    // Step 2: split per chromosome
    split_by_chr = split_by_subpop(preprocessed) |
		split_vcf_by_chr

	final_collection = merge_and_dedup(split_by_chr)
	
	pop_files = concat_chromosomes(final_collection)
	
	gwas_file = concat_subpops(pop_files)
}

process bcftools_filter {
    conda '/home/nicolas.lara/.conda/envs/imp_2'
    publishDir params.output_dir, mode: 'copy'

    input:
    path vcf

    output:
    tuple path("${sample}_filt.vcf.gz"), path("filtering_marker_table.txt")
	
    script:
    sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    bcftools view -i 'FORMAT/DP > ${params.depth} && 'FORMAT/GQ' >= ${params.quality}' ${vcf} -Oz -o DP_MAF_MISS_filter.vcf.gz
    bcftools view -m2 -M2 -v snps DP_MAF_MISS_filter.vcf.gz -Oz -o biallelic.vcf.gz
    bcftools index -c biallelic.vcf.gz
    bcftools view -t "^UNKNOWN" biallelic.vcf.gz -Oz -o "${sample}_filt.vcf.gz"

    ##create table of marker numbers
    echo -e "VCF File\tTotal Markers" > filtering_marker_table.txt
    for vcf in *.vcf.gz; do
        total_markers=\$(bcftools stats "\$vcf" | grep "^SN" | grep "number of records" | awk '{print \$6}')
        vcf_basename=\$(basename "\$vcf")
        echo -e "\${vcf_basename}\t\${total_markers}" >> filtering_marker_table.txt
    done
    """
}

process split_by_subpop {
    publishDir params.output_dir, mode: 'copy'
    conda '/home/nicolas.lara/.conda/envs/imp_2'

    input:
    path vcf         // Input VCF file
    path pop_table   // Population table file

    output:
    path "*.vcf.gz"

    script:
    """
    bgzip ${vcf}
    # Extract sample names from the VCF file header
    bcftools query -l ${vcf}.gz > all_samples.txt


    ##set internal field separator to comma
    IFS=","
    # For each population, subset out the population samples and parents
    while read Cross_ID Parent_1 Parent_2; do
        grep "^\${Cross_ID}-" all_samples.txt > "\${Cross_ID}_list.txt"
        echo "\${Parent_1}" >> "\${Cross_ID}_list.txt"
        echo "\${Parent_2}" >> "\${Cross_ID}_list.txt"
        bcftools view -S \${Cross_ID}_list.txt ${vcf}.gz -Oz -o \${Cross_ID}_subset.vcf.gz
    done < ${pop_table}
    """
}

process merge_and_dedup {
    input:
    tuple val(chrom), path(files)

    output:
    path("merged_${chrom}.tsv")
	
	script:
	"""
	library(data.table)
	library(glue)
	library(dplyr)
	
	##read in chrom VCFs as table
	m1 = data.frame(fread(cmd=glue("gzip -dc data/raw_vcf/{geno_list[1]} | grep -v '^##'")))
	m2 = data.frame(fread(cmd=glue("gzip -dc data/raw_vcf/{geno_list[2]} | grep -v '^##'")))
	m1[setdiff(colnames(m2), colnames(m1))] <- NA
	m2[setdiff(colnames(m1), colnames(m2))] <- NA
	m3 <- rbind(m1, m2[colnames(m1)])
	
	dup_samp <- table(m3$ID)
	dup_samp <- names(dup_samp[dup_samp>1])

	m3_orig <- m3[!(m3$ID %in% dup_samp),]

	get_consensus_call <- function(sample) {
	  calls <- tstrsplit(grep("\\./\\.", sample[[1]], invert = TRUE, value = TRUE), ":")
	  if (length(calls) == 0) {return(sample[[1]][1])}
	  GT <- names(sort(table(calls[1]), decreasing=T))[1]
	  GT_pos <- grep(GT, calls[[1]])
	  AD <- names(sort(table(calls[[2]][GT_pos]), decreasing=T))[1]
	  DP <- max(as.numeric(calls[[3]][GT_pos]))
	  GQ <- max(as.numeric(calls[[4]][GT_pos]))
	  PL <- names(sort(table(calls[[5]][GT_pos]), decreasing=T))[1]
	  return(paste(GT, AD, DP, GQ, PL, sep=":"))
	}

	dedup <- do.call(rbind, lapply(dup_samp, function(marker) {
	  mark_table <- m3[m3$ID == marker,]
	  ref <- names(sort(table(mark_table$REF), decreasing=T))[1]
	  rev_table <- mark_table[mark_table$REF != ref,]
	  rev_table[] <- lapply(rev_table, function(x) (gsub("0/", "n/", x)))
	  rev_table[] <- lapply(rev_table, function(x) (gsub("/0", "/n", x)))
	  rev_table[] <- lapply(rev_table, function(x) (gsub("1/", "0/", x)))
	  rev_table[] <- lapply(rev_table, function(x) (gsub("/1", "/0", x)))
	  rev_table[] <- lapply(rev_table, function(x) (gsub("n", "1", x)))
	  #swap ref and alt
	  temp <- rev_table$REF
	  rev_table$REF <- rev_table$ALT
	  rev_table$ALT <- temp
	  ##combine corrected row with majority ref rows
	  mark_table <- rbind(mark_table[mark_table$REF == ref,], rev_table)
	  concensus <- data.frame(c(mark_table[1, 1:9], apply(mark_table[10:ncol(mark_table)], 2, get_consensus_call)))
	  return(concensus)
	}))

	m3_new <- rbind(m3_orig, dedup) %>% arrange(X.CHROM, POS) %>% rename(CHROM = X.CHROM)


	##read in VCF
	m1_header <- readLines(paste0("data/raw_vcf/", geno_list[1]), n =100)
	m1_header <- m1_header[1:(grep("#CHROM", m1_header)-1)]
	m3_body_header <- paste0("#", paste(colnames(m3_new), collapse="\t"))
	vcf_output <- c(m1_header, m3_body_header, m3_new)
	writeLines(vcf_output, "${sample}_dedup.vcf")	
	"""
}

process beagle_impute {
    conda '/home/nicolas.lara/.conda/envs/imp_2'
    publishDir "${params.output_dir}", mode: 'copy'
    memory '100 GB'
    time '24h'


    input:
    path vcf

    output: 
    path "${sample}_imp.vcf.gz"

    script:
\    sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    bgzip ${vcf}
    beagle gt=${vcf}.gz \
            out=${sample}_imp \
            map=/project/guedira_seq_map/nico/SunRILs_population_description/linkage_map/monotonic/fam_GBS_monotonic.map \
            nthreads=40 \
            window=350
    """
}

process concat_subpops {
//not written yet
	input:
    path(tsv_files)

    output:
    path("${params.output_dir}/final_merged.tsv")

    script:
    """
    cat \$(ls merged_*.tsv | sort) > ${params.output_dir}/final_merged.tsv
    """
}