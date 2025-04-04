## run GWAS using ASRgwas package
## Nicolas A. H. Lara

library(here)
library(gaston)
library(dplyr)
library(GAPIT)
library(popkin)

setwd(here())


blues <- read.delim("data/blues.csv", sep=",") %>%
  rename(Taxa = Entry) %>%
  select(-c(Cross_ID)) %>%
  data.frame()

genotype <- read.vcf("data/SunRILs_prod_filt_imp.vcf.gz", convert.chr=F)
##Mapping vector for renaming chromosomes
mapping <- c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D",
             "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D",
             "7A", "7B", "7D")
genotype@snps$chr <- match(genotype@snps$chr, mapping)
genotype <- select.inds(genotype, id %in% blues$Taxa)

# kinship <- popkin(t(as.matrix(genotype)))
geno <- data.frame(Taxa = genotype@ped$id,
                   as.matrix(genotype)) %>%
  arrange(factor(Taxa, levels=blues$Taxa))


geno_map <- select(genotype@snps, c(id, chr, pos)) %>%
    rename(SNP = id, Chromosome = chr, Position = pos)

blues <- filter(blues, Taxa %in% geno$Taxa)

bonf_threshold <- (0.1 / ncol(genotype))

## run MLM model
for (trait in c("flowering", "Height", "Powdery_mildew", "WDR")) {
  print(trait)
  blu <- blues[,c("Taxa", trait)]
  gapit <- GAPIT(Y = blu, GD = geno, GM = geno_map, model= "MLM", file.output=F)
  gtable <- filter(gapit$GWAS, P.value <= bonf_threshold)
  gtable$trait <- trait
  gtable$model <- "MLM"
  if (exists("MLM_table")) {MLM_table  <- rbind(MLM_table , gtable)} else {MLM_table  <- gtable}
}
write.table(MLM_table, "outputs/gapit_mlm_gwas.tsv", quote=F, sep="\t", row.names=F)

## run BLINK model
for (trait in c("flowering", "Height", "Powdery_mildew", "WDR")) {
  blu <- blues[,c("Taxa", trait)]
  gapit <- GAPIT(Y = blu, GD = geno, GM = geno_map, model= "BLINK", file.output=F)
  gtable <- filter(gapit$GWAS, P.value <= bonf_threshold)
  gtable$trait <- trait
  gtable$model <- "BLINK"
  if (exists("BLINK_table")) {BLINK_table  <- rbind(BLINK_table , gtable)} else {BLINK_table  <- gtable}
}
write.table(BLINK_table, "outputs/gapit_blink_gwas.tsv", quote=F, sep="\t", row.names=F)


## run Farm-CPU model
for (trait in c("flowering", "Height", "Powdery_mildew", "WDR")) {
  blu <- blues[,c("Taxa", trait)]
  gapit <- GAPIT(Y = blu, GD = geno, GM = geno_map, model= "FarmCPU", file.output=F)
  gtable <- filter(gapit$GWAS, P.value <= bonf_threshold)
  gtable$trait <- trait
  gtable$model <- "FarmCPU"
  if (exists("FarmCPU_table")) {FarmCPU_table  <- rbind(FarmCPU_table , gtable)} else {FarmCPU_table  <- gtable}
}
write.table(FarmCPU_table, "outputs/gapit_farmcpu_gwas.tsv", quote=F, sep="\t", row.names=F)
