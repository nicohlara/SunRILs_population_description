##step 1: basic bcftools filtering.
##longterm goal: offload most of this step to this part
#bcftools view -i "FILTER/DP>1 && MAF>0.0012 && F_MISSING<0.1' in.vcf > out.vcf

library(gaston)

#setwd("C:/Users/nalara/Documents/GitHub/playground/filter_test")
setwd("/90daydata/guedira_seq_map/nico/SunRILs_adhoc")

vcf_file_name<-"SunRILs_prod_filt.vcf.gz"
geno <- gaston::read.vcf(vcf_file_name, convert.chr=F)

##filter SNP heterozygocity
geno <- select.snps(geno, hz <= 0.4)
##filter by excessive missing data per line
geno <- select.inds(geno, 
                    NAs/dim(geno)[2] < 0.7)

##rename lines
geno@ped$id <- sub("\\:.*", "", geno@ped$id)
geno@ped$famid <- ifelse(grepl("UX", geno@ped$id), sub("\\-.*", "", geno@ped$id), "Parent")

# geno_pop <- select.inds(geno, famid != "Parent")
# geno_par <- select.inds(geno, famid == "Parent")


# Function to calculate the most common value, its count, and NA count
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



if (identical(rownames(vcf_content), vcf_df$ID)) {
  vcf_content <- apply(vcf_content, 2, as.character)
  vcf_df <- cbind(vcf_df, vcf_content)
  vcf_body <- apply(vcf_df, 1, paste, collapse = "\t")
} else {
  print("Mismatch, check marker order")
}

vcf_body_head <- paste0("#", paste(colnames(vcf_df), collapse="\t"))
vcf_output <- c(vcf_head, vcf_body_head, vcf_body)
writeLines(vcf_output, "SunRILs_prod_filt2.vcf")
