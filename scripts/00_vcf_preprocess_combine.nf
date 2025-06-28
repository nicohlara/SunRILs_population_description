nextflow.enable.dsl=2

params.basedir = "project/guedira_seq_map/nico/SunRILs_population_description"
params.vcf_files = "${params.basedir}/data/raw_vcf"
params.population_table = "${params.basedir}/data/cross_info.csv"
params.output_dir = "${params.basedir}/data/processed_vcf"
params.monotonic_maps = "${params.basedir}/linkage_map/monotonic"
// bcftools filtering parameters
params.depth = '3'
params.quality = '20'

Channel
    .from(splitted_vcfs)
    .map { file -> 
        def cross_id = file.getBaseName().split('_')[1]
        tuple(cross_id, file)
    }
    .groupTuple()
	
workflow {
    // Collect input VCFs
    Channel.fromPath("${params.vcf_files}/*.vcf.gz")
        .set { raw_vcfs }

    // Broadcast the population table for pairing
    Channel.value(file(params.population_table))
        .set { pop_table_ch }

    // Step 1: Filter and split
    raw_vcfs
        | bcftools_filter
        | combine(pop_table_ch)
        | split_by_subpop
        | groupTupleBySubpop()
        | merge_and_dedup
        | beagle_impute
}

process bcftools_filter {
    input:
    path vcf

    output:
    tuple path("${sample}_filt.vcf.gz")
	
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
    input:
    tuple val(cross_id), path(vcf_files)

    output:
	tuple val(cross_id), path("${sample}_${cross_id}_filt.vcf.gz")
	
    script:
    sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    # Extract sample names from the VCF file header
    bcftools query -l ${vcf} > all_samples.txt

    ##set internal field separator to comma
    IFS=","
    # For each population, subset out the population samples and parents
    while read Cross_ID Parent_1 Parent_2; do
        grep "^\${Cross_ID}-" all_samples.txt > "\${Cross_ID}_list.txt"
        echo "\${Parent_1}" >> "\${Cross_ID}_list.txt"
        echo "\${Parent_2}" >> "\${Cross_ID}_list.txt"
        bcftools view -S \${Cross_ID}_list.txt ${vcf}.gz -Oz -o \${sample}_${Cross_ID}_subset.vcf.gz
    done < ${pop_table}
    """
}

process merge_and_dedup {
    input:
	tuple val(cross_id), path(vcf_files)
    
	output:
    tuple val(cross_id), path("${cross_id}_dedup.vcf")
	
	script:
	"""
	#!/usr/bin/env Rscript
	library(data.table)
	library(glue)
	library(dplyr)
	
	##read in subpop VCFs as table
	m1 = data.frame(fread(cmd=glue("gzip -dc data/raw_vcf/{vcf1} | grep -v '^##'")))
	m2 = data.frame(fread(cmd=glue("gzip -dc data/raw_vcf/{vcf2} | grep -v '^##'")))
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
	writeLines(vcf_output, "${cross_id}_dedup.vcf")	
	"""
}

process beagle_impute {
    publishDir "${params.output_dir}", mode: 'copy'
    memory '100 GB'
    time '24h'

    input:
	tuple val(cross_id), path(vcf)

    output: 
    path "${cross_id}_imp.vcf.gz"

    script:
    """
    bgzip ${vcf}
    beagle gt=${vcf}.gz \
            out=${cross_id}_imp \
            map=${params.monotonic}/${cross_id}_GBS_monotonic.map \
            nthreads=40 \
            window=350
    """
}