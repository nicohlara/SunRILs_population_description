nextflow.enable.dsl=2

//to run must have nextflow and miniconda3 loaded.
//see yaml of conda env packages included in project
// takes bcf and gaston filtered vcf files, combines them, then splits them into biparentals and imputes using beagle5 and biparental monotonic maps

//  Set up parameters
// tools
params.conda = '/home/nicolas.lara/.conda/envs/imputation'
// directories
params.vcf_dir = '/90daydata/guedira_seq_map/nico/SunRILs_filt_vcf'
params.output_dir = '/90daydata/guedira_seq_map/nico/SunRILs_filt_vcf/20250624_output'
params.pop_table = '/project/guedira_seq_map/nico/SunRILs_population_description/data/cross_info.csv'
params.linkage_maps = '/project/guedira_seq_map

// create output directory
new File(params.output_dir).mkdirs()



    



process gaston_clean {
    conda '/home/nicolas.lara/.conda/envs/imp_2'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple path(vcf), path(table)

    output:
    tuple path("${sample}_filt2.vcf"), path("filtering_marker_table.txt")
	
    script:
    sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    #!/usr/bin/env Rscript
    library(gaston)
    geno <- gaston::read.vcf("${vcf}", convert.chr=F)
	
    geno <- select.snps(geno, hz <= ${params.hz})

    ##filter by excessive missing data per line
    geno <- select.inds(geno, NAs/dim(geno)[2] < ${params.missing_line})

    comp_table <- read.delim("${table}")
    comp_table <- rbind(comp_table, c("gaston_filter", dim(geno)[2]))

	
    ##rename lines
    geno@ped\$id <- sub("\\\\:.*", "", geno@ped\$id)
    geno@ped\$famid <- ifelse(grepl("UX", geno@ped\$id), sub("\\\\-.*", "", geno@ped\$id), "Parent")
    tags <- c("-SMTISSUE", "-NWG", "-A+", "-A-", "+")
    geno@ped\$id <- gsub(paste0("\\\\b(", paste0(tags, collapse="|"), ")\\\\b"), "", geno@ped\$id)

    ##UX1999 and UX2031 are the same family with switched parents
    geno@ped\$id <- gsub("UX1999", "UX2031-99", geno@ped\$id)
	
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
    print("deduplicating samples")
    dup_samp <- table(geno@ped\$id)
    dup_samp <- names(dup_samp[dup_samp>1])
    for (ind in dup_samp) {
      print(ind)
      par <- as.matrix(select.inds(geno, id == ind))
      cs <- as.data.frame(t((sapply(par,consensus_info))))
      cs1 <- ifelse(cs\$count > nrow(par)/2, cs\$most_common, NA)
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
    print("Create new header")
    ##copy over header and add info
    vcf_head <- readLines("${vcf}", n=500)
    head_end <- grep("#CHROM", vcf_head) -1
    vcf_head <- vcf_head[1:head_end]

    vcf_df <- data.frame(CHROM=geno@snps\$chr,
						 POS=as.character(geno@snps\$pos),
						 ID=geno@snps\$id,
						 REF=geno@snps\$A1,
						 ALT=geno@snps\$A2,
						 QUAL=gsub("0", ".", geno@snps\$quality),
						 FILTER="PASS",
						 INFO=paste0("QualityScore=", geno@snps\$quality),
						 FORMAT="GT")

    vcf_content <- t(as.matrix(select.inds(geno, !(id %in% dup_samp))))
	if (exists("dedup")) {
      if (identical(colnames(dedup), row.names(vcf_content))) {
		vcf_content <- cbind(vcf_content, t(dedup))
	  } else {
	    print("Mismatch: check marker order")
      }
	}
	vcf_content[vcf_content==0] <- "0/0"
    vcf_content[vcf_content==1] <- "0/1"
    vcf_content[vcf_content==2] <- "1/1"
    vcf_content[is.na(vcf_content)] <- "./."
	
	
    #if (exists("dedup")) {
    #  if (identical(colnames(dedup), row.names(vcf_content))) {
  	#  vcf_content <- cbind(vcf_content, t(dedup))
	#  vcf_content[vcf_content==0] <- "0/0"
	#  vcf_content[vcf_content==1] <- "0/1"
	#  vcf_content[vcf_content==2] <- "1/1"
	#  vcf_content[is.na(vcf_content)] <- "./."
    #  } else {
	#  print("Mismatch: check marker order")
    #  }
    #}  else {
    #  vcf_content[vcf_content==0] <- "0/0"
    #  vcf_content[vcf_content==1] <- "0/1"
    #  vcf_content[vcf_content==2] <- "1/1"
    #  vcf_content[is.na(vcf_content)] <- "./."
    #}

    comp_table <- rbind(comp_table, c("gaston_consensus", dim(vcf_content)[1]))
	
    if (identical(rownames(vcf_content), vcf_df\$ID)) {
	  vcf_content <- apply(vcf_content, 2, as.character)
	  vcf_df <- cbind(vcf_df, vcf_content)
	  vcf_body <- apply(vcf_df, 1, paste, collapse = "\t")
    } else {
	  print("Mismatch, check marker order")
    }

    vcf_body_head <- paste0("#", paste(colnames(vcf_df), collapse="\t"))
    vcf_output <- c(vcf_head, vcf_body_head, vcf_body)
    writeLines(vcf_output, "${sample}_filt2.vcf")
    write.table(comp_table, file="filtering_marker_table.txt", quote=F, sep="\t", row.names=F)
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
    sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    bgzip ${vcf}
    beagle gt=${vcf}.gz \
            out=${sample}_imp \
            map=/90daydata/guedira_seq_map/nico/SunRILs/HPC-GBS-Pipeline/SynOp_RIL906_v1.0_GBS_monotonic.map \
            nthreads=40 \
            window=350
    """
}

workflow {
    Channel.fromPath(params.vcf_dir + '/*.vcf.gz')
        .set { vcf_files }
    Channel.fromPath(params.pop_table)
	.set { pop_table }


    vcf_raw = bcftools_filter(vcf_files) | 
	gaston_clean

    imp_vcf = vcf_raw.map { it[0] }
    beagle_impute(imp_vcf)
}

