nextflow.enable.dsl=2

//to run must have nextflow and miniconda3 loaded.
//see yaml of conda env packages included in project

//  Set up parameters
// tools
params.conda = '/home/nicolas.lara/.conda/envs/imputation'
// directories
//params.vcf_dir = '/90daydata/guedira_seq_map/nico/SunFilt'
params.vcf_dir = '/90daydata/guedira_seq_map/nico/SunFilt_round2'
params.output_dir = '/90daydata/guedira_seq_map/nico/SunFilt_round2/20250611_filter'
params.pop_table = '/project/guedira_seq_map/nico/SunRILs_population_description/data/cross_info.csv'

// BCFtools filtering parameters
params.depth = '3'
params.quality = '20'
params.MAF = '0.01'
params.missing = '0.4'

// Gaston filtering parameters
params.hz = '0.2'
params.missing_line = '0.5'

// create output directory
new File(params.output_dir).mkdirs()

process bcftools_filter {
    conda '/home/nicolas.lara/.conda/envs/imp_2'
    publishDir params.output_dir, mode: 'copy'

    input:
    path vcf

    output:
    tuple path("${sample}_filt.vcf.gz"), path("filtering_marker_table.txt")
	
    script:
//    sample = vcf.baseName
    sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    bcftools view -i 'FORMAT/DP > ${params.depth} && 'FORMAT/GQ' >= ${params.quality} && MAF > ${params.MAF} && F_MISSING < ${params.missing}' ${vcf} -Oz -o DP_MAF_MISS_filter.vcf.gz
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

process gaston_clean {
    // consider splitting this into 2 steps, an aggregation and output, then a gaston clean and output
    conda '/home/nicolas.lara/.conda/envs/imp_2'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple path(vcf), path(table)

    output:
    tuple path("${sample}_filt2.vcf"), path("filtering_marker_table.txt")
	
    script:
//    sample = vcf.baseName
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
      cs <- as.data.frame(t((apply(par, 2, consensus_info))))
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
	  vcf_content[vcf_content==0] <- "0/0"
	  vcf_content[vcf_content==1] <- "0/1"
	  vcf_content[vcf_content==2] <- "1/1"
	  vcf_content[is.na(vcf_content)] <- "./."
      } else {
	  print("Mismatch: check marker order")
      }
    }  else {
      vcf_content[vcf_content==0] <- "0/0"
      vcf_content[vcf_content==1] <- "0/1"
      vcf_content[vcf_content==2] <- "1/1"
      vcf_content[is.na(vcf_content)] <- "./."
    }

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
//    sample = vcf.baseName
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

process extract_samples {
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

workflow {
    Channel.fromPath(params.vcf_dir + '/*.vcf.gz')
        .set { vcf_files }
    Channel.fromPath(params.pop_table)
	.set { pop_table }


    vcf_raw = bcftools_filter(vcf_files) | 
	gaston_clean

    imp_vcf = vcf_raw.map { it[0] }
    beagle_impute(imp_vcf)


    // Extract sample names and create sample lists
//    sample_lists = extract_samples(imp_vcf, pop_table)

//    subset_vcf(sample_lists, vcf_raw)

    // Subset VCF based on sample lists
//    sample_lists.into { sample_files }
 //   sample_files.map { sample_list ->
  //      subset_vcf(sample_list, vcf_file)
   // }
}

