library(gaston)
library(chromoMap)
library(here)
library(dplyr)

setwd(here())

genotype <- read.vcf("data/SunRILs_prod_filt_imp.vcf.gz", convert.chr=F)


chr_file <- data.frame(chr = unique(genotype@snps$chr), start = 1)
a <- genotype@snps %>% group_by(chr) %>% summarise(max = max(pos))
chr_file <- merge(chr_file, a, by='chr')

anno_file <- genotype@snps %>% select(id, chr, pos) %>% mutate(end = pos+1)
