## merge vcfs, filter, and trim
## Nicolas A. H. Lara

##in bash:
# find -name "*_imp.vcf.gz" > list.txt
# while read file; do
#   sorted_file="${file/input_file/sorted_input_file}"
# 
#   bcftools index "$sorted_file"
#   bcftools sort "$file" -Oz -o "$sorted_file"
# 
#   # Extract marker IDs
#   base=$(basename "$file" .vcf.gz)
#   bcftools query -f '%ID\n' "$sorted_file" > "${base}_ids.txt"
# done < list.txt

# id_files <- grep("ids", list.files("data/processed_vcf/", full.names = T), value=T)
# id_lists <- lapply(id_files, function(file) {
#   read.delim(file, header = FALSE, col.names = "ID")$ID
# })
# id_counts <- unlist(id_lists) %>%
#   table() %>%
#   as.data.frame()
# colnames(id_counts) <- c("ID", "Count")
# id_majority <- id_counts %>%
#   filter(Count > 2*length(id_files) / 3)
# nrow(id_majority)
# write.table(id_majority$ID, "data/processed_vcf/well_covered_markers.txt", sep="\t", quote = F, row.names=F, col.names=F)
# 
# while read file; do
#   base=$(basename "$file" .vcf.gz)
#   bcftools view -i 'ID=@consensus_markers.txt' -Oz -o "${base}_prep.vcf.gz" "$file"
#   bcftools index "${base}_prep.vcf.gz"
# done < list.txt
# 
# bcftools merge *_prep.vcf.gz --force-samples -Oz -o SunRILs_imp.vcf.gz
# beagle gt=SunRILs_imp.vcf.gz out=SunRILs_mergeimp map=../../linkage_map/monotonic/consensus_GBS_monotonic.map nthreads=40 window=350

library(here)
library(gaston)
library(dplyr)
library(stringr)

setwd(here())

blues <- read.delim("data/blues.csv", sep=",")
genotype <- read.vcf("data/processed_vcf/SunRILs_imp_filt.vcf.gz", convert.chr=F)
pedigree <- read.csv("data/cross_info.csv", header=F, col.names = c("Cross_ID", "Parent_1", "Parent_2"))
genotype@ped$famid <- ifelse(grepl("UX", genotype@ped$id), str_sub(genotype@ped$id, 1, 6), 'Parent')


##merge parents
consensus_info <- function(column) {
  non_na_values <- column[!is.na(column)]
  if (length(non_na_values) == 0) {
    return(NA)
  }
  most_common <- as.numeric(names(sort(table(non_na_values), decreasing = TRUE)[1]))
  certainty <- sum(non_na_values == most_common)/(length(column)-length(column[is.na(column)]))
  if (certainty > 0.5) {
    return(most_common)
  } else {
    return(NA)
  }
}

for (parent in unique(c(pedigree$Parent_1, pedigree$Parent_2))) {
  parset <- grep(parent, genotype@ped$id, value=T)
  print(parset)
  if (length(parset) > 1) {
    cs <- matrix((apply(as.matrix(select.inds(genotype, id %in% parset)), 2, consensus_info)), nrow=1)
    colnames(cs) <- genotype@snps$id; rownames(cs) <- parent
    cs.fam <- data.frame(famid=parent,
                         id=parent,
                         father=0,
                         mother=0,
                         sex=0,
                         pheno=NA)
    cs.bim <- genotype@snps %>% select(chr, id, dist, pos, A1, A2)
    cs.bed <- as.bed.matrix(cs, cs.fam, cs.bim)
    geno <- select.inds(genotype, id %in% grep(parent, genotype@ped$id, value=T, invert=T))
    genotype <- gaston::rbind(geno, cs.bed)
  }
}

# ##select only lines with phenotypic data
genotype1 <- select.inds(genotype, id %in% blues$Entry)
# ##filter SNPs
genotype2 <- select.snps(genotype1, callrate == 1)
genotype3 <- select.snps(genotype2, N1/nrow(genotype2) < .5)
#filter lines
genotype4 <- select.inds(genotype3, N1/ncol(genotype3) < 0.5)
genotype5 <- select.inds(genotype4, NAs/ncol(genotype4) == 0)
genotype5 <- select.inds(genotype5, !(id %in% paste0("UX1992-", c(319, 8, 336, 273, 303, 9, 56, 211, 316 ))))
##LD thin
# genotype6 <- LD.thin(genotype5, threshold=0.8, max.dist=350e6)



dim(genotype); dim(genotype1); dim(genotype2); dim(genotype3); dim(genotype4); dim(genotype5); dim(genotype6)

write.bed.matrix(genotype5, "data/SunRILs_imp_filtmerge")
