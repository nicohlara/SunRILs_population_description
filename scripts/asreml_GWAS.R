## run GWAS using ASRgwas package
## Nicolas A. H. Lara

set.wd(here())

library(ASRgwas)
library(vcfR)

blues <- read.delim("data/blues.csv")
genotype <- read.vcfR("data/SunRILs_filt_imp.vcf.gz")

##filter genotype by entries in blues file

##filter blues by entries in genotype


##preprocess for gwas
gwas_obj <- pre.gwas(pheno.data = blues,
                     indiv = 'Genotype',
                     geno.data = geno.asr,
                     resp = 'flowering',
                     map.data = geno.map)




gwas <- gwas.asreml(pheno.data = gwas_obj$pheno.data,
                         resp = 'flowering', gen = 'Genotype',
                         Kinv = gwas_obj$Kinv,
                         geno.data = gwas_obj$geno.data, map.data = gwas_obj$map.data,
                         pvalue.thr = 0.0005, bonferroni = FALSE)

sign.markers <- gwas_test$gwas.sel$marker
geno.data.sel <- gwas_obj$geno.data[, sign.markers]
set <- select.marker(gwas.object = gwas_test,
                     geno.data.sel = geno.data.sel,
                     ref.vc = 1, pvalue.thr = 0.01)






pheno_subset <- filter(phenotype, Year %in% c(2022, 2023)) %>%
  rename(Genotype = Entry)
# pheno_subset <- filter(pheno_subset, Cross_ID == 'UX2013')
#library(ASRgenomics)
#pheno_subset <- rename(phenotype, Genotype = Entry)
# pheno_subset <- filter(phenotype, !is.na(height_range)) %>%
#   rename(Genotype = Entry, meas_height = Height, Height = height_range) %>%
#   mutate(Height = Height*100)
genotype_subset <- select.inds(genotype, id %in% pheno_subset$Genotype) 
pheno_subset <- filter(pheno_subset, Genotype %in% genotype_subset@ped$id) %>%
  select(c(Location, Year, row, column, Tray, Cross_ID, Genotype, flowering, Height, WDR, Powdery_mildew))
geno.asr <- as.matrix(genotype_subset)
geno.map <- genotype_subset@snps %>%
  dplyr::select(id, chr, pos) %>%
  rename(marker = id, chrom = chr)


library(ASRgwas)

gwas_obj <- pre.gwas(pheno.data = pheno_subset,
                     indiv = 'Genotype',
                     geno.data = geno.asr,
                     resp = 'flowering',
                     map.data = geno.map)
##tested using best naive guesses as to model factors
# gwas_test <- gwas.asreml(pheno.data = gwas_obj$pheno.data,
#             resp = 'flowering',
#             gen = 'Genotype',
#             Kinv = gwas_obj$Kinv,
#             fixedf = c('Location', 'Year'),
#             geno.data = gwas_obj$geno.data,
#             map.data = gwas_obj$map.data)\
gwas_test <- gwas.asreml(pheno.data = gwas_obj$pheno.data,
                         resp = 'flowering', gen = 'Genotype',
                         fixedf = c('Location', 'Year'), # residual = c('row', 'column'),
                         Kinv = gwas_obj$Kinv,
                         geno.data = gwas_obj$geno.data, map.data = gwas_obj$map.data,
                         pvalue.thr = 0.0005, bonferroni = FALSE)

sign.markers <- gwas_test$gwas.sel$marker
geno.data.sel <- gwas_obj$geno.data[, sign.markers]
set <- select.marker(gwas.object = gwas_test,
                     geno.data.sel = geno.data.sel,
                     ref.vc = 1, pvalue.thr = 0.01)
# set$gwas.sel
# set$call
# 
# qq.plot(gwas_test$gwas.all)

plot <- manhattan.plot(gwas.table=gwas_test$gwas.all, pvalue.thr = 1e-4, point.size = 3) +     theme(plot.background = element_rect(fill='transparent', color=NA))
plot
#ggsave(filename = paste0(dir, '../figures/', 'SunRILs_GWAS', '.png'), plot, width = 24, height = 8, dpi = 300, bg = 'transparent') 
pm_filter <- filter(gwas_test$gwas.sel, (effect >= .35 | effect <=-.35) & std.error < 0.1 & p.value < 1e-05)

ASRgwas::map.plot(map.data = gwas_test$gwas.all,
                  tag.markers = pm_filter$marker) + #gwas_test$gwas.sel$marker) +
  ggtitle('Powdery_mildew QTL')