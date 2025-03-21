nextflow.enable.dsl=2

//to run must have nextflow and miniconda3 loaded.
//see yaml of conda env packages included in project

//  Set up parameters
// tools
params.conda = '/home/nicolas.lara/.conda/envs/imputation'
// directories
params.vcf_dir = '/90daydata/guedira_seq_map/nico/SunFilt'
params.output_dir = '/90daydata/guedira_seq_map/nico/SunFilt'

// BCFtools filtering parameters
params.depth = '1'
params.MAF = '0.02'
params.missing = '0.2'

// Gaston filtering parameters
params.hz = .2
params.missing_line = 0.6


process bcftools_filter {
//    conda '${params.conda}'
    conda '/home/nicolas.lara/.conda/envs/imputation'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path vcf

    output:
    tuple path("${sample}_filt.vcf.gz"), path("filtering_marker_table.txt"), emit "${sample}_filt.vcf.gz"
	
    script:
    sample = vcf.baseName
    """
    bcftools view -i 'FORMAT/DP>${params.depth} && MAF > ${params.MAF} && F_MISSING<${params.missing}' ${vcf} -Oz -o DP_MAF_MISS_filter.vcf.gz
    bcftools view -m2 -M2 -v snps DP_MAF_MISS_filter.vcf.gz -Oz -o  biallelic.vcf.gz
    bcftools view -t "^UNKNOWN" biallelic.vcf.gz -Oz -o  "${sample}_filt.vcf.gz"
    echo -e "VCF File\tTotal Markers" > filtering_marker_table.txt
    for v in *.vcf.gz; do
        total_markers=\$(bcftools stats "\$v" | grep "^SN" | grep "number of records" | awk '{print \$6}')
        vcf_basename=\$(basename "\$v")
        echo -e "\${vcf_basename}\t\${total_markers}" >> filtering_marker_table.txt
    done
    """
}

process gaston_clean {
	// consider splitting this into 2 steps, an aggregation and output, then a gaston clean and output
    conda '/home/nicolas.lara/.conda/envs/imputation'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple path(vcf), path(table)

    output:
    path "${sample}_filt2.vcf.gz"
	
    script:
    sample = vcf.baseName
    """
	library(gaston)
	geno <- gaston::read.vcf(${vcf}, convert.chr=F)
	
	geno <- select.snps(geno, hz <= ${params.hz})
    ##filter by excessive missing data per line
    geno <- select.inds(geno, NAs/dim(geno)[2] < ${params.missing_line})
	
	comp_table <- read.delim(${table})
	comp_table <- rbind(comp_table, c("gaston_filter", dim(geno)[2]))
	
	##rename lines
	geno@ped$id <- sub("\\:.*", "", geno@ped$id)
	geno@ped$famid <- ifelse(grepl("UX", geno@ped$id), sub("\\-.*", "", geno@ped$id), "Parent")
	consensus_info <- function(column) {
	  non_na_values <- column[!is.na(column)]
	  if (length(non_na_values) == 0) {
		return(c(most_common = NA, count = 0))
	  }
	  most_common <- as.numeric(names(sort(table(non_na_values), decreasing = TRUE)[1]))
	  return(c(most_common = most_common, 
			   count = sum(non_na_values == most_common)))
	}

	##get duplicated samples
	dup_samp <- table(geno@ped$id)
	dup_samp <- names(dup_samp[dup_samp>1])
	for (ind in dup_samp) {
	  print(ind)
	  par <- as.matrix(select.inds(geno, id == ind))
	  cs <- as.data.frame(t((apply(par, 2, consensus_info))))
	  cs1 <- ifelse(cs$count > nrow(par)/2, cs$most_common, NA)
	  if (exists("dedup")) { 
		cs2 <- data.frame(t(cs1))
		row.names(cs2) <- ind; colnames(cs2) <- colnames(par)
		dedup <- rbind(dedup, cs2)
	  } else {
		  names(cs1) <- row.names(cs)
		  dedup <- data.frame(t(cs1))
		  row.names(dedup) <- ind
	  }
	}
	###create new VCF
	##copy over header and add info
	vcf_head <- readLines(vcf_file_name, n=500)
	head_end <- grep("#CHROM", vcf_head) -1
	vcf_head <- vcf_head[1:head_end]
	# vcf_hqlesead[head_end+1] <- "##R/gaston select.snps hz <= 0.05, select.inds missing < 0.7"

	vcf_df <- data.frame(CHROM=geno@snps$chr,
						 POS=as.character(geno@snps$pos),
						 ID=geno@snps$id,
						 REF=geno@snps$A1,
						 ALT=geno@snps$A2,
						 QUAL=gsub("0", ".", geno@snps$quality),
						 FILTER=geno@snps$chr,
						 INFO=paste0("QualityScore=", geno@snps$quality),
						 FORMAT="GT")

	vcf_content <- t(as.matrix(select.inds(geno, !(id %in% dup_samp))))
	if (identical(colnames(dedup), row.names(vcf_content))) {
	  vcf_content <- cbind(vcf_content, t(dedup))
	  vcf_content[vcf_content==0] <- "0/0"
	  vcf_content[vcf_content==1] <- "0/1"
	  vcf_content[vcf_content==2] <- "1/1"
	  vcf_content[is.na(vcf_content)] <- "./."
	} else {
	  print("Mismatch: check marker order")
	}

	comp_table <- rbind(comp_table, c("gaston_consensus", dim(vcf_content)[1]))
	
	if (identical(rownames(vcf_content), vcf_df$ID)) {
	  vcf_content <- apply(vcf_content, 2, as.character)
	  vcf_df <- cbind(vcf_df, vcf_content)
	  vcf_body <- apply(vcf_df, 1, paste, collapse = "\t")
	} else {
	  print("Mismatch, check marker order")
	}

	vcf_body_head <- paste0("#", paste(colnames(vcf_df), collapse="\t"))
	vcf_output <- c(vcf_head, vcf_body_head, vcf_body)
	writeLines(vcf_output, "${sample}_filt2.vcf.gz")
	write.table(comp_table, file="filtering_marker_table.txt", quote=F, sep="\t", row.names=F)
	"""
}

workflow {
    Channel.fromPath(params.vcf_dir + '/*.vcf.gz') \
        .set { vcf_files }

	bcftools_filter(vcf_files) |
		gaston_clean
}
