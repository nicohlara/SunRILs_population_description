## run GWAS using ASRgwas package
## Nicolas A. H. Lara

library(here)
library(ASRgwas)
library(gaston)
library(dplyr)

setwd(here())


blues <- read.delim("data/blues.csv", sep=",") %>%
  mutate(Cross_ID = as.factor(Cross_ID)) %>%
  rename(Genotype = Entry)
genotype <- read.vcf("data/SunRILs_prod_filt_imp.vcf.gz", convert.chr=F)

##take smallest biparental and get proportion of whole pop
##then multiply by reasonable MAF for F4 genotyped individuals
table(blues$Cross_ID)
MAF_threshold <- (110/sum(table(blues$Cross_ID)))*0.5

##filter genotype by entries in blues file
# genotype <- select.inds(genotype, id %in% blues$Entry) 
# Mapping vector for renaming chromosomes
# mapping <- c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", 
#              "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", 
#              "7A", "7B", "7D")
# genotype@snps$chr <- match(genotype@snps$chr, mapping)

  
##filter blues by entries in genotype
# blues <- filter(blues, Entry %in% genotype@ped$id) %>%
#   rename(Genotype = Entry)

geno.map <- genotype@snps %>%
  dplyr::select(id, chr, pos) %>%
  rename(marker = id, chrom = chr)

bonf_threshold <- (0.1 / ncol(genotype))

for (trait in c("flowering", "Height", "Powdery_mildew", "WDR")) {
  ##preprocess for gwas
  # gb <- blues[,c("Genotype", trait)]
  gwas_obj <- pre.gwas(pheno.data = blues,
                       indiv = 'Genotype',
                       geno.data = as.matrix(genotype),
                       resp = trait,
                       map.data = geno.map,
                       maf=MAF_threshold)
  
  gwas <- gwas.asreml(pheno.data = gwas_obj$pheno.data,
                      resp = trait,
                      gen = 'Genotype',
                      Kinv = gwas_obj$Kinv,
                      geno.data = gwas_obj$geno.data,
                      map.data = gwas_obj$map.data, 
                      bonferroni = F,
                      pvalue.thr=bonf_threshold)
  
  gtable <- gwas$gwas.sel
  gtable$trait <- trait
  if (exists("GWAS_table")) {GWAS_table <- rbind(GWAS_table, gtable)} else {GWAS_table <- gtable}
}

write.table(GWAS_table, "outputs/asreml_gwas.tsv", quote=F, sep="\t", row.names=F)
