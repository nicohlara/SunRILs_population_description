## merge vcfs, filter, and trim
## Nicolas A. H. Lara

library(here)
library(gaston)
library(dplyr)
library(stringr)

setwd(here())

blues <- read.delim("data/blues.csv", sep=",")
genotype <- read.vcf("data/processed_vcf_allsequencedata_20250817/SunRILs_imp_filt.vcf.gz", convert.chr=F)
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
    cs.fam <- data.frame(famid="Parent",
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
hist(genotype2@snps$N1/nrow(genotype2))
genotype3 <- select.snps(genotype2, N1/nrow(genotype2) < .2) #could go to .2 without issue
#filter lines
hist(genotype3@ped$N1/ncol(genotype3))
genotype4 <- select.inds(genotype3, N1/ncol(genotype3) < 0.15) #could go to 0.1 or 0.15 without issue
genotype5 <- select.inds(genotype4, NAs/ncol(genotype4) == 0)
genotype5 <- select.inds(genotype5, !(id %in% paste0("UX1992-", c(319, 8, 336, 273, 303, 9, 56, 211, 316 ))))
##LD thin
# genotype6 <- LD.thin(genotype5, threshold=0.8, max.dist=350e6)



dim(genotype); dim(genotype1); dim(genotype2); dim(genotype3); dim(genotype4); dim(genotype5)

write.bed.matrix(genotype2, "data/SunRILs_imp_filtmerge")



g <- read.vcf("data/processed_vcf/SunRILs_raw.vcf.gz", convert.chr=F)
g@ped$famid <- ifelse(grepl("UX", g@ped$id), str_sub(g@ped$id, 1, 6), 'Parent')

# g1 <- select.snps(g, maf > 0.011)
gt <- select.inds(g, famid == 'Parent')
gp <- select.snps(gt, callrate > 0.8)
g1 <- select.snps(g, id %in% gp@snps$id)
g1 <- select.snps(g1, N1/nrow(g1) < 0.1) 
g1 <- select.inds(g1, N1/ncol(g1) < 0.1)

variation_marker <- c()
for (fam in grep("UX", unique(g1@ped$famid), value=T)) {
  print(fam)
  gt <- select.inds(g1, famid == fam)
  gt <- select.snps(gt, NAs/nrow(gt) < 0.5) 
  g2 <- select.snps(g1, id %in% gt@snps$id)
  gt <- select.snps(gt, maf > 0.3)
  variation_marker <- unique(c(variation_marker, gt@snps$id))
}
g3 <- select.snps(g2, id %in% variation_marker)
print(dim(g2)); print(dim(g3))

depth <- read.delim('depth.txt', header=F)
depth1 <- as.data.frame(lapply(depth[-c(1:3)], function(x) sub("^.*?:.*?:(.*?):.*?:.*$", "\\1", x)))
dv <- as.numeric(as.vector(as.matrix(depth1)))
hist(dv[dv < 100])
