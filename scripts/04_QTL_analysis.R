
library(here)
library(dplyr)
library(tidyr)
library(glue)
library(data.table)
library(gaston)
library(ggplot2)
library(ggnewscale)
setwd(here())

pedigree <- read.csv("data/cross_info.csv", header=F, col.names = c("Cross_ID", "Parent_1", "Parent_2"))

# 
# 
# ## HEATMAP OF MARKER EFFECTS
# pks <- read.delim("outputs/qtl2_cim.tsv") %>%
#   filter(trait == 'PM') %>%
#   mutate(phys_pos = round(as.numeric(sapply(strsplit(Pos, "_"), "[", 2))/1e6, 1),
#          chr = gsub("\\..", "", chr)) %>%
#   rename(chromosome = chr) %>%
#   arrange(chromosome, phys_pos) %>%
#   rename(Cross_ID = cross)
# 
# ##divide into S/L
# cls <-  read.delim("data/chrom_lengths_split.txt")
# pks$chromosome <- apply(pks, 1, function(x) {y <- as.numeric(x[[17]])*1e6; cls[cls$Chromosome == x[[6]] & cls$Start <= y
#                                                                                & cls$End > y,'Arm']})
# 
# ##add in GWAS summary
# gwas_results <- read.delim("outputs/combined_GWAS.tsv") %>%
#   mutate(lod = -log10(P.value))
# gwas_results$chromosome <- apply(gwas_results, 1, function(x) {y <- as.numeric(x[[2]]); cls[cls$Chromosome == x[[1]] & cls$Start <= y
#                                                                                             & cls$End > y,'Arm']})
# gwas_results <- select(gwas_results, c(Position, Chromosome, chromosome, trait, lod)) %>%
#   filter(trait == 'PM') %>%
#   mutate(Cross_ID = 'GWAS', marker = glue("S{Chromosome}_{Position}") )
# genotype <- read.vcf("data/SunRILs_prod_filt_imp.vcf.gz", convert.chr=F)
# genotype@ped$family <- ifelse(grepl("UX", genotype@ped$id), str_sub(genotype@ped$id, 1, 6), 'Parent')
# genotype <- select.snps(genotype, id %in% gwas_results$marker)
# blues <- read.delim("data/blues.csv", sep=",") %>%
#   mutate(Cross_ID = as.factor(Cross_ID)) %>% 
#   rename(id = Entry)
# genotype <- select.inds(genotype, id %in% blues$id)
# 
# geno <- cbind(id = genotype@ped$id, as.matrix(genotype)) %>%
#   data.frame() 
# geno_long <-data.frame(geno) %>%
#   inner_join(blues[,c('id', 'PM')], by = "id") %>%
#   pivot_longer(cols = -c(id, PM), names_to = "marker", values_to = "genotype")
# 
# # 2. Compute means per marker/genotype
# geno_summary <- geno_long %>%
#   group_by(marker, genotype) %>%
#   summarise(mean_blue = mean(PM, na.rm = TRUE), .groups = "drop") %>%
#   pivot_wider(names_from = genotype, values_from = mean_blue, names_prefix = "mean_pheno_")
# 
# # View the result
# geno_centered <- geno_summary %>%
#   mutate(
#     midpoint = (mean_pheno_0 + mean_pheno_2) / 2,
#     adj_pheno_0 = mean_pheno_0 - midpoint,
#     adj_pheno_2 = mean_pheno_2 - midpoint
#   ) %>%
#   select(marker, adj_pheno_0, adj_pheno_2)
# 
# gwas <- merge(gwas_results, geno_centered, by= 'marker', all=T) %>%
#   rename(Pos = marker, AA = adj_pheno_0, BB = adj_pheno_2)
# 
# pks1 <- merge(pks, gwas, by=intersect(names(pks), names(gwas)), all=T)
# 
# peak_summary <- as.data.table(pks1)[, {
#   idx <- which.max(abs(AA))
#   .(
#     effect_peak = pos[idx],
#     phys_pos = phys_pos[idx],
#     LOD = lod[idx],
#     AA = AA[idx],
#     BB = BB[idx],
#     NC08 = NC08[idx],
#     HILLIARD = HILLIARD[idx],
#     GA13LE6 = GA13LE6[idx],
#     nmar = uniqueN(pos)
#   )
# }, by = .(trait, chromosome, Cross_ID)]
# 
# peak_summary1 <- merge(peak_summary, pedigree, by='Cross_ID') %>%
#   mutate(cross_label = paste(Cross_ID, Parent_2, Parent_1, sep=" - "))
# # peak_label = paste0(chromosome, ":", peak_start, "-", peak_end)) 
# 
# #convert to NAM founder allele coding
# peak_summary1 <- peak_summary1 %>%
#   rename("NC08-23383" = NC08,
#          "GA06493-13LE6" = GA13LE6) %>%
#   rowwise() %>%
#   mutate(
#     genotype_val = get(Parent_2),
#     missing_geno = ifelse(is.na(genotype_val), 1, 0),
#     # AA = if (!is.na(genotype_val) && genotype_val == 2) -1 * AA else AA,
#     # BB = if (!is.na(genotype_val) && genotype_val == 2) -1 * BB else BB
#   ) %>%
#   ungroup() %>%
#   select(-genotype_val)
# 
# peak_summary2 <- filter(peak_summary, Cross_ID == 'GWAS') %>%
#   mutate(cross_label = 'GWAS')
# peak_summary <- merge(peak_summary2, peak_summary1, 
#                       by=intersect(names(peak_summary2), names(peak_summary1)), all=T) %>%
#   mutate(chromosome = factor(chromosome, levels = cls$Arm))
# 
# 
# 
# ##cleaning up missing values manually
# # peak_summary[peak_summary$chromosome == '7AL' & peak_summary$Cross_ID == "UX1992", c("AA", "BB", "missing_geno")] <- 
# #   c(peak_summary[peak_summary$chromosome == '7AL' & peak_summary$Cross_ID == "UX1992", c("AA", "BB")]*-1, 0)
# 
# parentage_order <- pedigree[order(pedigree$Parent_2, decreasing=F),] %>% 
#   mutate(cross_label = paste(Cross_ID,Parent_2, Parent_1, sep=" - ")) %>%
#   select(cross_label) %>% as.vector
# peak_summary <- mutate(peak_summary, cross_label = factor(cross_label, levels = c( "GWAS", rev(parentage_order[[1]]))))
# 
# ggplot(peak_summary, aes(x = chromosome, y = cross_label, fill = AA)) +
#   geom_point(aes(fill = AA, size = LOD), shape = 22, colour = 'gray90') +
#   
#   # Create a factor for the legend and map it to shape
#   geom_point(
#     data = subset(peak_summary, missing_geno == 1),
#     aes(shape = "Missing parent genotype"),
#     size = 2
#   ) +
#   
#   scale_fill_gradient2(
#     low = "blue", mid = "white", high = "red", 
#     midpoint = 0, na.value = "grey60", name = "AA Effect"
#   ) +
#   scale_size(range = c(1, 10), name = 'LOD') +
#   scale_shape_manual(
#     name = NULL,
#     values = c("Missing parent genotype" = 21)
#   ) +
#   
#   theme_minimal(base_size = 12) +
#   labs(x = "Chromosome", y = "Cross") +
#   theme(
#     axis.text.x = element_text(angle = 60, hjust = 1),
#     axis.text.y = element_text(hjust = 0),
#     legend.box = "vertical",
#     legend.position = "right"
#   )
# 
# ggsave(file="figures/QTL_heatmap_chromArm.png", width=12, height=5)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ### CREATE QTL HEATMAP WITH EFFECTS
# pks <- read.delim("outputs/qtl2_cim.tsv") %>%
#   mutate(phys_pos = round(as.numeric(sapply(strsplit(Pos, "_"), "[", 2))/1e6, 1),
#          chr = gsub("\\..", "", chr)) %>%
#   rename(chromosome = chr) %>%
#   arrange(chromosome, phys_pos) %>%
#   rename(Cross_ID = cross)
# cls <-  read.delim("data/chrom_lengths_split.txt")
# pks$chromosome <- apply(pks, 1, function(x) {y <- as.numeric(x[[17]])*1e6; cls[cls$Chromosome == x[[6]] & cls$Start <= y
#                                                                                & cls$End > y,'Arm']})
# ##read in GWAS summary
# gwas_results <- read.delim("outputs/combined_GWAS.tsv") %>%
#   mutate(lod = -log10(P.value))
# gwas_results$chromosome <- apply(gwas_results, 1, function(x) {y <- as.numeric(x[[2]]); cls[cls$Chromosome == x[[1]] & cls$Start <= y
#                                                                                             & cls$End > y,'Arm']})
# gwas_results <- select(gwas_results, c(Position, Chromosome, chromosome, trait, lod)) %>%
#   mutate(Cross_ID = 'GWAS', marker = glue("S{Chromosome}_{Position}") )
# genotype <- read.bed.matrix("data/SunRILs_imp_filtmerge")
# genotype <- select.snps(genotype, id %in% gwas_results$marker)
# blues <- read.delim("data/blues.csv", sep=",") %>%
#   mutate(Cross_ID = as.factor(Cross_ID)) %>% 
#   rename(id = Entry)
# genotype <- select.inds(genotype, id %in% blues$id)
# geno <- cbind(id = genotype@ped$id, as.matrix(genotype)) %>%
#   data.frame() 
# 
# plot_qtl_heatmap <- function(pks, gwas_results, geno, blues, trt) {
#   selection <- c('id', trt)
#   geno_long <-data.frame(geno) %>%
#     inner_join(blues[,selection], by = "id") %>%
#     pivot_longer(cols = -c(selection), names_to = "marker", values_to = "genotype")
#   
#   # 2. Compute means per marker/genotype
#   geno_summary <- geno_long %>%
#     group_by(marker, genotype) %>%
#     dplyr::summarise(mean_blue := mean(.data[[trt]], na.rm = TRUE), .groups = "drop") %>%
#     pivot_wider(names_from = genotype, values_from = mean_blue, names_prefix = "mean_pheno_")
#   
#   # View the result
#   geno_centered <- geno_summary %>%
#     mutate(
#       midpoint = (mean_pheno_0 + mean_pheno_2) / 2,
#       adj_pheno_0 = mean_pheno_0 - midpoint,
#       adj_pheno_2 = mean_pheno_2 - midpoint
#     ) %>%
#     select(marker, adj_pheno_0, adj_pheno_2)
#   
#   gwas <- merge(dplyr::filter(gwas_results, trait == trt), geno_centered, by= 'marker', all=T) %>%
#     rename(Pos = marker, AA = adj_pheno_0, BB = adj_pheno_2)
#   
#   pks1 <- merge(pks, gwas, by=intersect(names(pks), names(gwas)), all=T) %>%
#     filter(trait == trt)
#   
#   peak_summary <- as.data.table(pks1)[, {
#     idx <- which.max(abs(AA))
#     .(
#       effect_peak = pos[idx],
#       phys_pos = phys_pos[idx],
#       LOD = lod[idx],
#       AA = AA[idx],
#       BB = BB[idx],
#       NC08 = NC08[idx],
#       HILLIARD = HILLIARD[idx],
#       GA13LE6 = GA13LE6[idx],
#       nmar = uniqueN(pos)
#     )
#   }, by = .(trait, chromosome, Cross_ID)]
#   
#   peak_summary1 <- merge(peak_summary, pedigree, by='Cross_ID') %>%
#     mutate(cross_label = paste(Cross_ID, Parent_2, Parent_1, sep=" - "))
#   
#   #convert to NAM founder allele coding
#   peak_summary1 <- peak_summary1 %>%
#     rename("NC08-23383" = NC08,
#            "GA06493-13LE6" = GA13LE6) %>%
#     rowwise() %>%
#     mutate(
#       genotype_val = get(Parent_2),
#       missing_geno = ifelse(is.na(genotype_val), 1, 0),
#     ) %>%
#     ungroup() %>%
#     select(-genotype_val)
#   
#   peak_summary2 <- filter(peak_summary, Cross_ID == 'GWAS') %>%
#     mutate(cross_label = 'GWAS')
#   peak_summary <- merge(peak_summary2, peak_summary1, 
#                         by=intersect(names(peak_summary2), names(peak_summary1)), all=T) %>%
#     mutate(chromosome = factor(chromosome, levels = cls$Arm))
#   
#   parentage_order <- pedigree[order(pedigree$Parent_2, decreasing=F),] %>% 
#     mutate(cross_label = paste(Cross_ID,Parent_2, Parent_1, sep=" - ")) %>%
#     select(cross_label) %>% as.vector
#   peak_summary <- mutate(peak_summary, cross_label = factor(cross_label, levels = c( "GWAS", rev(parentage_order[[1]]))))
#   
#   ggplot(peak_summary3, aes(x = chromosome, y = cross_label, fill = AA)) +
#     geom_point(aes(fill = AA, size = LOD), shape = 22, colour = 'gray90') +
#     
#     # Create a factor for the legend and map it to shape
#     geom_point(
#       data = subset(peak_summary3, missing_geno == 1),
#       aes(shape = "Missing parent genotype"),
#       size = 2
#     ) +
#     
#     scale_fill_gradient2(
#       low = "blue", mid = "white", high = "red", 
#       midpoint = 0, na.value = "grey60", name = "AA Effect"
#     ) +
#     scale_size(range = c(1, 8), name = 'LOD') +
#     scale_shape_manual(
#       name = NULL,
#       values = c("Missing parent genotype" = 21)
#     ) +
#     
#     theme_minimal(base_size = 12) +
#     labs(x = "Chromosome", y = "Cross", title = glue("QTL effect for {trt}")) +
#     theme(
#       axis.text.x = element_text(angle = 60, hjust = 1),
#       axis.text.y = element_text(hjust = 0),
#       legend.box = "vertical",
#       legend.position = "right"
#     )
#   
#   ggsave(file=glue("figures/QTL_heatmap_chromArm_{trt}.png"), width=12, height=5)
# }
# 
# for (trt in c('WDR', "HD", "PM", "Height")) {
#   plot_qtl_heatmap(pks, gwas_results, geno, blues, trt)
# }
# 
# 
# 
# 

### QTL EFFECT HEATMAP
##read in data and initialization
cls <-  read.delim("data/chrom_lengths_split.txt")
qtl2_cim <- read.delim("outputs/qtl2_cim.tsv") %>%
  mutate(Position = round(as.numeric(sapply(strsplit(Pos, "_"), "[", 2))/1e6, 1),
         chr = gsub("\\..", "", chr)) %>%
  rename(Cross_ID = cross, marker = Pos, Chromosome = chr) %>%
  select(marker, Cross_ID, trait, Chromosome, Position, lod, AA, BB, intercept, NC08, HILLIARD, GA13LE6)
gwas_results <- read.delim("outputs/combined_GWAS.tsv") %>%
  mutate(lod = -log10(P.value),
         Cross_ID = 'GWAS',
         marker = glue("S{Chromosome}_{Position}"),
         Position = round(Position/1e6, 1)) %>%
  select(-c(P.value, model))
blues <- read.delim("data/blues.csv", sep=",") %>%
  rename(id = Entry)
genotype <- read.bed.matrix("data/SunRILs_imp_filtmerge")
geno <- select.snps(genotype, id %in% gwas_results$marker)
geno <- select.inds(geno, id %in% blues$id)
geno <- cbind(id = geno@ped$id, as.matrix(geno)) %>%
  data.frame() 

##massage into correct structures
geno_long <- geno %>%
  inner_join(blues[,-1], by = "id") %>%
  pivot_longer(cols = -colnames(blues)[-1], names_to = "marker", values_to = "genotype")

for (trt in colnames(blues[-c(1:2)])) {
  geno_summary <- geno_long %>%
    group_by(marker, genotype) %>%
    dplyr::summarise(mean_blue := mean(.data[[trt]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = genotype, values_from = mean_blue, names_prefix = "mean_pheno_")
  
  # View the result
  geno_centered <- geno_summary %>%
    mutate(
      midpoint = (mean_pheno_0 + mean_pheno_2) / 2,
      adj_pheno_0 = mean_pheno_0 - midpoint,
      adj_pheno_2 = mean_pheno_2 - midpoint,
      trait = trt) %>%
    select(marker, adj_pheno_0, adj_pheno_2, trait) %>%
    filter(!is.na(adj_pheno_0))
  if (exists("gece")) {gece <- rbind(gece, geno_centered)} else {gece <- geno_centered}
}

gwas <- merge(gwas_results, gece, by = c("marker", "trait")) %>%
  rename(AA = adj_pheno_0, BB = adj_pheno_2)

pks1 <- merge(qtl2_cim, gwas, by=intersect(names(qtl2_cim), names(gwas)), all=T)
pks1$Chromosome <- apply(pks1, 1, function(x) {y <- as.numeric(x[[5]])*1e6; cls[cls$Chromosome == x[[4]] & cls$Start <= y
                                                                                & cls$End > y,'Arm']})
                                                                               

peak_summary <- as.data.table(pks1)[, {
  idx <- which.max(abs(AA))
  .(
    Position = Position[idx],
    LOD = lod[idx],
    AA = AA[idx],
    BB = BB[idx],
    NC08 = NC08[idx],
    HILLIARD = HILLIARD[idx],
    GA13LE6 = GA13LE6[idx],
    nmar = uniqueN(marker)
  )
}, by = .(trait, Chromosome, Cross_ID)]

peak_summary1 <- merge(peak_summary, pedigree, by='Cross_ID') %>%
  mutate(cross_label = paste(Cross_ID, Parent_2, Parent_1, sep=" - "))

#convert to NAM founder allele coding
peak_summary1 <- peak_summary1 %>%
  rename("NC08-23383" = NC08,
         "GA06493-13LE6" = GA13LE6) %>%
  rowwise() %>%
  mutate(
    genotype_val = get(Parent_2),
    missing_geno = ifelse(is.na(genotype_val), 1, 0),
  ) %>%
  ungroup() %>%
  select(-genotype_val)

peak_summary2 <- filter(peak_summary, Cross_ID == 'GWAS') %>%
  mutate(cross_label = 'GWAS')

peak_summary3 <- merge(peak_summary2, peak_summary1, 
                      by=intersect(names(peak_summary2), names(peak_summary1)), all=T) %>%
  mutate(Chromosome = factor(Chromosome, levels = cls$Arm))

parentage_order <- pedigree[order(pedigree$Parent_2, decreasing=F),] %>% 
  mutate(cross_label = paste(Cross_ID,Parent_2, Parent_1, sep=" - ")) %>%
  select(cross_label) %>% as.vector
peak_summary3 <- mutate(peak_summary3, cross_label = factor(cross_label, levels = c( "GWAS", rev(parentage_order[[1]]))))

plot_QTL_heatmap <- function(peaks, trt) {
  pks <- filter(peaks, trait == trt)
  ggplot(pks, aes(x = Chromosome, y = cross_label, fill = AA)) +
    geom_point(aes(fill = AA, size = LOD), shape = 22, colour = 'gray90') +
    
    # Create a factor for the legend and map it to shape
    geom_point(
      data = subset(pks, missing_geno == 1),
      aes(shape = "Missing parent genotype"),
      size = 2
    ) +
    
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", 
      midpoint = 0, na.value = "grey60", name = "AA Effect"
    ) +
    scale_size(range = c(1, 8), name = 'LOD') +
    scale_shape_manual(
      name = NULL,
      values = c("Missing parent genotype" = 21)
    ) +
    
    theme_minimal(base_size = 12) +
    labs(x = "Chromosome", y = "Cross", title = glue("QTL effect for {trt}")) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      axis.text.y = element_text(hjust = 0),
      legend.box = "vertical",
      legend.position = "right"
    )
  ggsave(file=glue("figures/QTL_heatmap_chromArm_{trt}.png"), width=12, height=5)
}

for (trt in c('WDR', "HD", "PM", "Height")) {
  plot_QTL_heatmap(peak_summary3, trt)
}

wdrhd <- peak_summary3 %>%
  select(Chromosome, cross_label, trait, AA, LOD) %>%
  filter(trait %in% c("WDR", "HD")) %>%
  mutate(
    x = as.numeric(droplevels(Chromosome)),
    y = as.numeric(as.factor(cross_label)),
    triangle = if_else(trait == "WDR", "topleft", "bottomright")
  )

make_triangle_coords <- function(x, y, triangle_type, lod, lod_range = c(1, 8)) {
  # Normalize LOD to 0–1 range
  lod_scaled <- (lod - lod_range[1]) / diff(lod_range)
  lod_scaled <- pmax(pmin(lod_scaled, 1), 0)  # clamp between 0 and 1
  
  # triangle size radius (range: 0.2 to 0.5)
  size <- 0.2 + lod_scaled * 0.3
  
  if (triangle_type == "topleft") {
    tibble(
      x1 = c(x - size, x - size, x + size),
      y1 = c(y + size, y - size, y + size)
    )
  } else {
    tibble(
      x1 = c(x + size, x + size, x - size),
      y1 = c(y - size, y + size, y - size)
    )
  }
}
lod_min <- min(wdrhd$LOD, na.rm = TRUE)
lod_max <- max(wdrhd$LOD, na.rm = TRUE)
triangles <- wdrhd %>%
  rowwise() %>%
  mutate(coords = list(make_triangle_coords(x, y, triangle, LOD, lod_range = c(lod_min, lod_max)))) %>%
  unnest(coords) %>%
  ungroup()

ggplot() +
  # Trait A: green-red
  geom_polygon(
    data = filter(triangles, trait == "WDR"),
    aes(x = x1, y = y1, group = interaction(Chromosome, cross_label, trait), fill = AA),
    color = "gray80"
  ) +
  scale_fill_gradient2(low = "#ff0055", mid = "white", high = "#5500ff", midpoint = 0, name = 'WDR') +
  new_scale_fill() +
  
  # Trait B: blue-orange
  geom_polygon(
    data = filter(triangles, trait == "HD"),
    aes(x = x1, y = y1, group = interaction(Chromosome, cross_label, trait), fill = AA),
    color = "gray80"
  ) +
  scale_fill_gradient2(low = "#0055ff", mid = "white", high = "#ff5500", midpoint = 0, name = 'HD') +
  
  scale_x_continuous(breaks = unique(triangles$x), labels = unique(wdrhd$Chromosome)) +
  scale_y_continuous(breaks = unique(triangles$y), labels = unique(wdrhd$cross_label)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text.y = element_text(hjust = 0)
  ) +
  labs(x = "Chromosome", y = "Cross", title = "Split QTL Effects")
ggsave(file=glue("figures/QTL_heatmap_WDR_HD_comb.png"), width=12, height=5)
