## merge vcfs, filter, and trim
## Nicolas A. H. Lara

##in bash:
# find -name "*_imp.vcf.gz" > list.txt
# cat list.txt | while read file; do
#   bcftools index "${file/input_file/sorted_input_file}"
#   bcftools sort $file -Oz -o "${file/input_file/sorted_input_file}"
# done
# bcftools merge *_imp.vcf.gz --force-samples -Oz -o SunRILs_imp.vcf.gz
# beagle gt=SunRILs_imp.vcf.gz out=SunRILs_mergeimp map=../../linkage_map/monotonic/consensus_GBS_monotonic.map nthreads=40 window=350

library(here)
library(gaston)
library(dplyr)

setwd(here())

blues <- read.delim("data/blues.csv", sep=",")
genotype <- read.vcf("data/processed_vcf/SunRILs_imp.vcf.gz", convert.chr=F)
pedigree <- read.csv("data/cross_info.csv", header=F, col.names = c("Cross_ID", "Parent_1", "Parent_2"))


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
    genotype <- rbind(geno, cs.bed)
  }
}

##select only lines with phenotypic data
genotype <- select.inds(genotype, id %in% blues$Entry)
##filter SNPs
genotype <- select.snps(genotype, callrate > 0.6)
genotype <- select.snps(genotype, N1/nrow(genotype) < .5)
#filter lines
genotype <- select.inds(genotype, N1/ncol(genotype) < 0.5)
genotype <- select.inds(genotype, NAs/ncol(genotype) < 0.5)

write.bed.matrix(genotype, "data/SunRILs_imp_filtmerge")
