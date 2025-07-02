nextflow.enable.dsl=2

params.basedir = "/project/guedira_seq_map/nico/SunRILs_population_description"
params.vcf_files = "${params.basedir}/data/raw_vcf"
params.population_table = "${params.basedir}/data/cross_info.csv"
params.output_dir = "${params.basedir}/data/processed_vcf"
params.monotonic = "${params.basedir}/linkage_map/monotonic"
// bcftools filtering parameters
params.depth = '3'
params.quality = '20'

workflow {
    new File(params.output_dir).mkdirs()

    // Load input VCFs
    Channel.fromPath("${params.vcf_files}/*.vcf.gz")
        .set { raw_vcfs }

    // Broadcast population table
    Channel.value(file(params.population_table))
        .set { pop_table_ch }

    // Step 1: Filter and split
    filtered_vcfs = bcftools_filter(raw_vcfs)

    // Step 2: Combine each filtered VCF with the population table
    vcf_with_pop_table = filtered_vcfs.combine(pop_table_ch)
        | prefilter_populations

    // Step 3: Split filtered VCFs by subpopulations
    subpop_vcfs = split_by_subpop(vcf_with_pop_table)
        .flatten()
        .map { file -> 
            def cross_id = file.getBaseName().tokenize("_")[-2]
            tuple(cross_id, file)
        }

    // Step 4: Group subset VCFs by Cross_ID
    grouped_subpop_vcfs = subpop_vcfs
        .map { cross_id, file -> tuple(cross_id, file) }
        .groupTuple()

    // Step 5: Split into paired and singleton Cross_IDs
    //paired_ch = grouped_subpop_vcfs
    //    .filter { cross_id, files -> files.size() == 2 }

    //singleton_ch = grouped_subpop_vcfs
    //    .filter { cross_id, files -> files.size() == 1 }

    // Step 6a: Merge and deduplicate
    //deduped_vcfs = paired_ch 
    deduped_vcfs = grouped_subpop_vcfs 
        .filter { cross_id, files -> files.size() in [1,2] }
        .map { cross_id, files -> tuple(cross_id, files) }
        | merge_and_dedup
        | compress_vcf

    // Step 6b: Pass singleton files directly
    //singleton_passthrough = singleton_ch
    //    .map { cross_id, files -> tuple(cross_id, files[0]) }

    // Step 7: Beagle imputation
    deduped_vcfs
      //  .mix(singleton_passthrough)
      //  .flatten()
        .map { cross_id, vcf ->
            println "Sending ${vcf} for imputation of ${cross_id}"
            tuple(cross_id, vcf)  // ensure tuple structure
        }
        | beagle_impute
}

process bcftools_filter {
    input:
    path vcf

    output:
    path "${sample}_filt.vcf.gz"
	
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

process clean_vcf_sample_ids {
    input:
    path vcf

    output:
    path "${sample}_rename.vcf.gz"
                                                                                                                        script:                                                                                                             sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    # Step 1: Extract original sample names
    bcftools query -l ${vcf} > original_samples.txt

    # Step 2: Clean sample names and prepare rename map
    paste original_samples.txt original_samples.txt | awk -v OFS='\\t' '
    {
        old = \$2

        # Step 1: Remove colon suffix
        sub(/:.*/, "", old)

        # Step 2: Remove unwanted tags
        gsub(/-SMTISSUE|-NOPLANT|-LOST|-NOTISSUE|-NWG|-A\\+|-A\\-|\\+/, "", old)

        # Step 3: Rename UX1999
        if (old == "UX1999") old = "UX2031-99"

        print \$1, old
    }' > rename_map.txt

    # Step 3: Reheader VCF using cleaned names
    bcftools reheader -s rename_map.txt -o ${sample}_rename.vcf.gz ${vcf}
    """
}

process prefilter_populations {
    input:
    tuple path(vcf), path(pop_table)

    output:
    tuple path(vcf), path("subset_pop_table.txt")

    script:
    """
    # Extract sample names from the VCF file header
    bcftools query -l ${vcf} > all_samples.txt

    # Step 2: Extract Cross_IDs from sample names
    grep '^UX[0-9]\\{4\\}-' all_samples.txt | sed 's/-.*//' | sort -u > cross_ids_present.txt

    # Step 3: Filter the population table to include only those Cross_IDs
    awk -F',' 'NR==FNR {keep[\$1]; next} \$1 in keep' cross_ids_present.txt ${pop_table} > subset_pop_table.txt
    """
}

process split_by_subpop {
    errorStrategy 'ignore'

    input:
    tuple path(vcf), path(pop_table)

    output:
    path("*_subset.vcf.gz")
	
    script:
    sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    # Extract sample names from the VCF file header
    bcftools query -l ${vcf} > all_samples.txt

    ##set internal field separator to comma
    IFS=","
    # For each population, subset out the population samples and parents
    while read Cross_ID Parent_1 Parent_2; do
        grep "^\${Cross_ID}-" all_samples.txt > temp_list.txt
        echo "\${Parent_1}" >> temp_list.txt
        echo "\${Parent_2}" >> temp_list.txt
        
        grep -Fxf temp_list.txt all_samples.txt > "\${Cross_ID}_list.txt"

        bcftools view -S \${Cross_ID}_list.txt ${vcf} -Oz -o ${sample}_\${Cross_ID}_subset.vcf.gz
    done < ${pop_table}

    ls -l *.vcf.gz || echo "No matching VCFs created."
    """
}

process merge_and_dedup {
    input:
    //tuple val(cross_id), path(vcf1), path(vcf2)
    tuple val(cross_id), path(vcfs)    

    output:
    tuple val(cross_id), path("${cross_id}_dedup.vcf")
	
	script:
	"""
	#!/usr/bin/env Rscript
	library(data.table)
	library(glue)
	library(dplyr)

	print("${cross_id}")

        tags <- c("-SMTISSUE", "-NWG", "-A+", "-A-", "+")
	
	##read in subpop VCFs as tables
	vcf_files <- c(${vcfs.collect { "'${it}'" }.join(", ")})
	m1 = data.frame(fread(cmd=glue("gzip -dc {vcf_files[1]} | grep -v '^##'")))
        colnames(m1) <- sub("\\\\:.*", "", colnames(m1))
        colnames(m1) <-  gsub(paste0("\\\\b(", paste0(tags, collapse="|"), ")\\\\b"), "", colnames(m1))
        colnames(m1) <- gsub("UX1999", "UX2031-99", colnames(m1))
	
	if (length(vcf_files) == 2) {
		m2 <- data.frame(fread(cmd = glue("gzip -dc {vcf_files[2]} | grep -v '^##'")))
	        colnames(m2) <- sub("\\\\:.*", "", colnames(m2))
	        colnames(m2) <-  gsub(paste0("\\\\b(", paste0(tags, collapse="|"), ")\\\\b"), "", colnames(m2))
	        colnames(m2) <- gsub("UX1999", "UX2031-99", colnames(m2))	

		m1[setdiff(colnames(m2), colnames(m1))] <- NA
		m2[setdiff(colnames(m1), colnames(m2))] <- NA
		m3 <- rbind(m1, m2[colnames(m1)])
	} else {
		m3 <- m1
	}
	
	print("Full dedup dimensions")
	print(dim(m3))

	dup_samp <- table(m3\$ID)
	dup_samp <- names(dup_samp[dup_samp>1])
	
	print("Duplicated samples")
	print(length(dup_samp))	

	m3_orig <- m3[!(m3\$ID %in% dup_samp),]

	get_consensus_call <- function(sample) {
	  sample <- sample[!is.na(sample)]
	  calls <- tstrsplit(grep("\\\\./\\\\.", sample, invert = TRUE, value = TRUE), ":")
	  if (length(calls) == 0) {return(sample[1])}
	  GT <- names(sort(table(calls[1]), decreasing=T))[1]
	  GT_pos <- grep(GT, calls[[1]])
	  AD <- names(sort(table(calls[[2]][GT_pos]), decreasing=T))[1]
	  DP <- max(as.numeric(calls[[3]][GT_pos]))
	  GQ <- max(as.numeric(calls[[4]][GT_pos]))
	  PL <- names(sort(table(calls[[5]][GT_pos]), decreasing=T))[1]
	  return(paste(GT, AD, DP, GQ, PL, sep=":"))
	}
	
	if (length(dup_samp) > 0) {
		dedup <- do.call(rbind, lapply(dup_samp, function(marker) {
		  mark_table <- m3[m3\$ID == marker,]
		  ref <- names(sort(table(mark_table\$REF), decreasing=T))[1]
		  rev_table <- mark_table[mark_table\$REF != ref,]
		  rev_table[] <- lapply(rev_table, function(x) (gsub("0/", "n/", x)))
		  rev_table[] <- lapply(rev_table, function(x) (gsub("/0", "/n", x)))
		  rev_table[] <- lapply(rev_table, function(x) (gsub("1/", "0/", x)))
		  rev_table[] <- lapply(rev_table, function(x) (gsub("/1", "/0", x)))
		  rev_table[] <- lapply(rev_table, function(x) (gsub("n", "1", x)))
		  #swap ref and alt
		  temp <- rev_table\$REF
		  rev_table\$REF <- rev_table\$ALT
		  rev_table\$ALT <- temp
		  ##combine corrected row with majority ref rows
		  mark_table <- rbind(mark_table[mark_table\$REF == ref,], rev_table)
		  print("Marker table dimensions")
		  print(dim(mark_table))
		  concensus <- data.frame(c(mark_table[1, 1:9], apply(mark_table[10:(ncol(mark_table)-1)], 2, get_consensus_call), mark_table[ncol(mark_table)]))
		  return(concensus)
		}))
	
		m3_new <- rbind(m3_orig, dedup)
	} else {
		m3_new <- m3_orig 
	}

	m3_new <- m3_new %>% arrange(X.CHROM, POS) %>% rename(CHROM = X.CHROM)
	print("Final vcf dimensions")
	print(dim(m3_new))
	##read in VCF
	m1_header <- readLines(vcf_files[1], n =100)
	m1_header <- m1_header[1:(grep("#CHROM", m1_header)-1)]
	m3_body_header <- paste0("#", paste(colnames(m3_new), collapse="\t"))
	m3_body <- apply(m3_new, 1, function(row) paste(trimws(as.character(row)), collapse = "\t"))
	print(class(m1_header))
        print(class(m3_body_header))
        print(class(m3_new))
	vcf_output <- c(m1_header, m3_body_header, m3_body)
	writeLines(vcf_output, "${cross_id}_dedup.vcf")	
	"""
}

process compress_vcf {
    input:
    tuple val(cross_id), path(vcf)

    output:
    tuple val(cross_id), path("*.gz")

    script:
    """
    bgzip -c ${vcf} > ${vcf.simpleName}.vcf.gz
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
    beagle gt=${vcf} \
            out=${cross_id}_imp \
            map=${params.monotonic}/${cross_id}_GBS_monotonic.map \
            nthreads=40 \
            window=350
    """
}
