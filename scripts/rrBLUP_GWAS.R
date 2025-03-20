## run GWAS using ASRgwas package
## Nicolas A. H. Lara

library(here)
library(rrBLUP)
library(gaston)
library(dplyr)

setwd(here())


blues <- read.delim("data/blues.csv", sep=",") %>%
  mutate(Cross_ID = as.factor(Cross_ID))
genotype <- read.vcf("data/SunRILs_prod_filt_imp.vcf.gz", convert.chr=F)

##take smallest biparental and get proportion of whole pop
##then multiply by reasonable MAF for F4 genotyped individuals
table(blues$Cross_ID)
MAF_threshold <- (110/sum(table(blues$Cross_ID)))*0.5^4 

##filter genotype by entries in blues file
genotype <- select.inds(genotype, id %in% blues$Entry)
# Mapping vector for renaming chromosomes
# mapping <- c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D",
#              "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D",
#              "7A", "7B", "7D")
# genotype@snps$chr <- match(genotype@snps$chr, mapping)
geno <- genotype@snps[,c("id", "chr", "pos")]
geno <- cbind(geno, t(as.matrix(genotype))-1)
##filter blues by entries in genotype
blues <- filter(blues, Entry %in% genotype@ped$id) %>%
  rename(Genotype = Entry)

# geno.map <- genotype@snps %>%
#   dplyr::select(id, chr, pos) %>%
#   rename(marker = id, chrom = chr)


bonf_threshold <- (0.1 / ncol(genotype))

for (trait in c("flowering", "Height", "Powdery_mildew", "WDR")) {
  ##preprocess for gwas
  pheno <- blues[,c("Genotype", trait)]
  gwas <- rrBLUP::GWAS(pheno=pheno,
                       geno=geno,
                       min.MAF=MAF_threshold,
                       #plot=F
                       )
  
  # gtable <- gwas$gwas.sel
  colnames(gwas)[4] <- "p.value"
  gwas$p.value <- 10^-gwas$p.value
  gwas <- filter(gwas, p.value <= bonf_threshold)
  gwas$trait <- trait
  if (exists("GWAS_table")) {GWAS_table <- rbind(GWAS_table, gwas)} else {GWAS_table <- gwas}
}

# GWAS_table$p.value <- 10^-GWAS_table$p.value
write.table(GWAS_table, "outputs/rrBLUP_gwas.tsv", quote=F, sep="\t", row.names=F)
