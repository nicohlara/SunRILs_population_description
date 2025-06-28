library(gaston)
library(glue)
geno_list <- c("SunRILs_production.vcf", "UX1992_reseq_production.vcf")

# m1 <- data.frame(as.matrix(read.vcf(file=paste0("data/", geno_list[1]), convert.chr = F)))
# m2 <- data.frame(as.matrix(read.vcf(file=paste0("data/", geno_list[2]), convert.chr = F)))
# 
# 
# ##read in VCF
# m1 <- readLines(paste0("data/", geno_list[1]))
# m1 <- m1[grep("#CHROM", m1):length(m1)]
# m1 <- gsub("#", "", m1[1])
# m2 <- fread(m1, sep="\t")

## read in vcf files
m1 = data.frame(fread(cmd=glue("gzip -dc data/raw_vcf/{geno_list[2]} | grep -v '^##'")), nrows=1000)#, sep='\t', header = TRUE, skip = '#CHROM', nrows=0))
m1 <- fread(cmd = glue("bcftools view -r 1A data/raw_vcf/{geno_list[1]} | grep -v '^##'"))
m1 <- fread(cmd = glue("awk '$1 == \"1A\" || $1 ~ /^#CHROM/' data/raw_vcf/{geno_list[1]}"))

cmd_str <- glue("awk '$1 == \"1A\" || $1 ~ /^#CHROM/' data/raw_vcf/{geno_list[1]}")
m1 <- fread(cmd = cmd_str)


lines <- readLines(glue("data/raw_vcf/{geno_list[1]}"))
lines_chr1A <- lines[grepl("^1A\\t|^#CHROM", lines)]
m1 <- fread(text = lines_chr1A)



m2 = data.frame(fread(cmd=glue("gzip -dc data/raw_vcf/{geno_list[2]} | grep -v '^##'")))#, sep='\t', header = TRUE, skip = '#CHROM', nrows=0))

## make vcf files equivalent and merge
m1[setdiff(colnames(m2), colnames(m1))] <- NA
m2[setdiff(colnames(m1), colnames(m2))] <- NA
m3 <- rbind(m1, m2[colnames(m1)])

## identify duplicated markers and merge



# a <- data.frame(matrix(data=c('1A', 1, 'T', 'A', 0, 0,
#                    '1A', 2, 'A', 'T', 0, 0,
#                    '1A', 3, 'T', 'A', 0, 0),
#             nrow=3, byrow = T,
#             dimnames=list(c(1,2,3), c('CHROM', 'POS', "REF", 'ALT', 'A', 'B'))))
# 
# b <- data.frame(matrix(data=c('1A', 1, 'T', 'A', 1, 1,
#                    '1A', 2, 'T', 'A', 1, 1,
#                    '1A', 4, 'T', 'A', 1, 1),
#             nrow=3, byrow = T,
#             dimnames=list(c(1,2,4), c('CHROM', 'POS', "REF", 'ALT', 'A', 'C'))))
# 
# a[setdiff(colnames(b), colnames(a))] <- NA
# b[setdiff(colnames(a), colnames(b))] <- NA
# d <- rbind(a, b[colnames(a)])

consensus_info <- function(column) {
  non_na_values <- column[!is.na(column)]
  if (length(non_na_values) == 0) {
    return(c(most_common = NA, count = 0))
  }
  most_common <- as.numeric(names(sort(table(non_na_values), decreasing = TRUE)[1]))
  return(c(most_common = most_common, 
           count = sum(non_na_values == most_common)))
}


consensus_call <- function(column) {
  calls <- column[column != "./." & !is.na(column)]
  if (length(unique(calls)) == 1) {return(calls)} else {return(NA)}
  ###CHANGE RETURN NA TO RETURN './.' AFTER TESTING
}


dup_samp <- table(m3$ID)
dup_samp <- names(dup_samp[dup_samp>1])
for (ind in dup_samp) {
  par <- m3[m3$ID == ind,]
  if (par[1, 'REF'] == par[2, 'REF'] & par[1, 'ALT'] == par[2, 'ALT']) {
    cs <- as.data.frame(t((sapply(par[10:ncol(par)],consensus_call))))
  } else {
    print(ind)
    print("REF/ALT mismatch")
  }
  if (exists("dedup")) {dedup <- rbind(dedup, cbind(par[1, 1:9], cs)) } else {dedup <- cbind(par[1, 1:9], cs)}
}







# 
# b <- m1[order(m1$ID) & m1$ID %in% a,] 
# d <- m2[order(m2$ID) & m2$ID %in% a,] 
# 
# 
# # g1 <- read.vcf(file=paste0("data/", geno_list[1]), convert.chr = F)
# # g2 <- read.vcf(file=paste0("data/", geno_list[2]), convert.chr = F)
# 
# m1[setdiff(colnames(m2), colnames(m1))] <- NA
# m2[setdiff(colnames(m1), colnames(m2))] <- NA
# # m1$id <- row.names(m1)
# # m2$id <- row.names(m2)
# 
# m3 <- rbind(m1, m2[colnames(m1)])























print("deduplicating samples")
consensus_info <- function(column) {
  non_na_values <- column[!is.na(column)]
  if (length(non_na_values) == 0) {
    return(c(most_common = NA, count = 0))
  }
  most_common <- as.numeric(names(sort(table(non_na_values), decreasing = TRUE)[1]))
  return(c(most_common = most_common, 
           count = sum(non_na_values == most_common)))
}

dup_samp <- table(m3$id)
dup_samp <- names(dup_samp[dup_samp>1])
for (ind in dup_samp) {
  print(ind)
  # par <- as.matrix(select.inds(geno, id == ind))
  par <- m3[m3$id == ind,]
  par <- par[!(colnames(par) == "id")]
  cs <- as.data.frame(t((sapply(par,consensus_info))))
  cs1 <- ifelse(cs$count >= nrow(par)/2, cs$most_common, NA)
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


vcf_content <- t(m3[!(m3$id %in% dup_samp),])
vcf_content <- vcf_content[!(row.names(vcf_content) == 'id'),]


# m3 <- m3[!(colnames(m3) == "id")]

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

print("Create new header")
##copy over header and add info
vcf_head <- readLines(paste0("data/", geno_list[1]), n=500)
head_end <- grep("#CHROM", vcf_head) -1
vcf_head <- vcf_head[1:head_end]

vcf_df <- data.frame(CHROM= sapply(strsplit(row.names(vcf_content), "_"), "[", 1),
                     POS=sapply(strsplit(row.names(vcf_content), "_"), "[", 2),
                     ID=row.names(vcf_content),
                     REF=geno@snps$A1,
                     ALT=geno@snps$A2,
                     QUAL=gsub("0", ".", geno@snps\$quality),
                     FILTER="PASS",
                     INFO=paste0("QualityScore=", geno@snps$quality),
                     FORMAT="GT")










###working mostly
library(dplyr)
library(glue)
m3

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
