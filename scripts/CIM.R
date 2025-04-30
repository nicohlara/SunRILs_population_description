## run CIM using r/qtl2
## Nicolas A. H. Lara

library(here)
library(qtl)
library(qtl2)
library(dplyr)

setwd(here())


pedigree <- read.csv("data/cross_info.csv")
blues <- read.delim("data/blues.csv", sep=",") %>%
  rename(genotype = Entry) %>%
  select(-Cross_ID)

# cross_id <- 'UX1992'
# chrom = '6A'
blues <- read.csv("2023_phenotype/trait_BLUEs.csv") %>%
  rename(genotype = Entry)
for (fam in pedigree$Cross_ID[-1]) {
  print(fam)
  dir <- glue("linkage_map/{fam}_map")
  cross <- read.cross(format='csv', file=glue("{dir}/{fam}_qtl2_map.csv"),
                      estimate.map=F, genotypes=c('AA', 'AB', 'BB'),
                      crosstype='riself')
  # if ('6A' %in% row.names(summary.map(cross))) {chrom = '6A'} else {chrom='6A.2'}
  cross <- subset(cross, ind = (cross$pheno$genotype %in% blues$genotype))
  cross$pheno <- merge(cross$pheno[1], blues, by='genotype', all=F)
  SunCross <- convert2cross2(cross)
  
  ##get map and calculate QTL probabilities
  map <- SunCross$gmap
  ##calculate QTL probabilities using geno prob
  pr <- calc_genoprob(SunCross, error_prob=0.002, cores=4)
  bonf_threshold <- -log10((0.1 / nrow(out)))
  ##calculate a kinship matrix of relationship among individuals
  kinship <- calc_kinship(pr, type='overall', use_allele_probs = FALSE)
  ##Perform genome scan using Haley-Knott regression (without kinship) or LMM (with kinship)
  for (trait in colnames(blues[-1])) {
    out <- scan1(pr, SunCross$pheno[, trait], kinship=kinship, cores=0)
    pks <- find_peaks(out, map, peakdrop=0.5, expand2markers=T, drop=0.5, threshold = 1.5)
    if (nrow(pks) > 0 ) {
      Pos <- c()
      for (i in 1:nrow(pks)) {
        Pos <- c(Pos, find_marker(map, pks[i,'chr'], pos=pks[i,'pos']))
      }
      pks$Pos <- Pos
      pks$cross <- fam; pks$trait <- trait
      if (exists('peaks')) {peaks <- rbind(peaks, pks)} else {peaks <- pks}
    }
    png(filename=glue("{dir}/{fam}_{trait}_CIM.png"),    
        width=750*3, height=300*3, res=72*3,
        bg='transparent')
    par(mar=c(4,4,6,1))
    plot(out, map, lodcolumn = 'pheno1', main=paste0("CIM of ", fam), bg='transparent')
    abline(h=bonf_threshold, col='#CC0000')
    dev.off()
  }
}
write.table(peaks, file="outputs/qtl2_cim.tsv", quote=F, sep="\t", row.names=F)

