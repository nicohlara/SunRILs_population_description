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
# blues <- read.csv("2023_phenotype/trait_BLUEs.csv") %>%
#   rename(genotype = Entry)
for (fam in pedigree$Cross_ID) {
  print(fam)
  # dir <- glue("linkage_map/{fam}_map")
  # cross <- read.cross(format='csv', file=glue("{dir}/{fam}_qtl2_map.csv"),
  #                     estimate.map=F, genotypes=c('AA', 'AB', 'BB'),
  #                     crosstype='riself')
  cross <- read.cross(format="csv",file=glue("linkage_map/maps/{fam}_linkage_map.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
  ##need to fix this in creation, not sure how it happened
  cross$pheno <- cross$pheno[,1:6]
  colnames(cross$pheno) <- gsub("\\.x", "", colnames(cross$pheno))
  # cross <- subset(cross, ind = (cross$pheno$genotype %in% blues$genotype))
  # cross$pheno <- merge(cross$pheno[1], blues, by='genotype', all=F)
  SunCross <- convert2cross2(cross)
  
  ##get map and calculate QTL probabilities
  map <- SunCross$gmap
  ##calculate QTL probabilities using geno prob
  pr <- calc_genoprob(SunCross, map = map, cores=4)
  bonf_threshold <- -log10((0.1 / totmar(cross)))
  ##calculate a kinship matrix of relationship among individuals
  # kinship <- calc_kinship(pr, type='overall', use_allele_probs = FALSE)
  ##Perform genome scan using Haley-Knott regression (without kinship) or LMM (with kinship)
  for (trait in colnames(blues[-1])) {
    out <- scan1(pr, SunCross$pheno[, trait])
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
    # png(filename=glue("{dir}/{fam}_{trait}_CIM.png"),   
    png(filename=glue("figures/CIM_plots/{fam}_{trait}_CIM.png"),
        width=750*3, height=300*3, res=72*3,
        bg='transparent')
    par(mar=c(4,4,6,1))
    plot(out, map, lodcolumn = 'pheno1', main=paste0("CIM of ", fam), bg='transparent')
    abline(h=bonf_threshold, col='#CC0000')
    dev.off()
  }
}
write.table(peaks, file="outputs/qtl2_cim.tsv", quote=F, sep="\t", row.names=F)


# 
# 
# ## qtl test
# cross <- jittermap(cross)
# cross <- calc.genoprob(cross, step=5)
# # cross$pheno <- cross$pheno[5]
# operm <- scanone(cross, method="hk", n.perm=1000)
# summary(operm, alpha = c(0.1, 0.2, 0.9, 1))
# lod_thresh <- summary(operm, alpha=0.95)[1]
# 
# 
# plot(scanone(cross, method="hk"))
# 
# qtl_result <- stepwiseqtl(cross, method="hk", additive.only = T, penalties=c(lod_thresh, 0, 0))
# plot(qtl_result)
# 
# cim_results <- cim(cross, n.marcovar=5, method="hk")
# plot(cim_results)
# 
# 
# ### CS1.1 test
# 
# # 
# # cross_old <- read.cross(format="csv",file="G:/My Drive/Nico_PhD.lnk/data/genotype/final_Gmap/UX1992n_GMap.csv",
# #                     estimate.map=FALSE, na.strings=c("-","NA"),
# #                     genotypes=c("A", "H", "B"), crosstype = "riself")
# #                     # genotypes=c("AA","H","BB"), crosstype="riself")
# # cross_old <- subset(cross_old, ind = (cross_old$pheno$genotype %in% blues$genotype))
# # cross_old$pheno <- merge(cross_old$pheno[1], blues, by='genotype', all=F)
# # # cross_old <- subset(cross_old, ind = (cross_old$pheno$genotype %in% oldblue$genotype))
# # # cross_old$pheno <- merge(cross_old$pheno[1], oldblue, by='genotype', all=F)
# # 
# # # 
# # # b <- merge(cross$pheno[,c("genotype", "Height")], cross_old$pheno[,c("genotype", "Height")], by="genotype")
# # # cor(b$Height.x, b$Height.y)
# # 
# # 
# # cross_old <- jittermap(cross_old)
# # cross_old <- calc.genoprob(cross_old, step=5)
# # plot(scanone(cross_old, pheno.col=4, method="hk"))
# # plot(scanone(cross, pheno.col=4, method="hk"))
# # 
# # # cross_old$pheno <- cross_old$pheno[5]
# # operm <- scanone(cross_old, pheno.col=5, method="hk", n.perm=1000)
# # summary(operm, alpha = c(0.1, 0.2, 0.9, 1))
# # lod_thresh <- summary(operm, alpha=0.5)[1]
# # # qtl_result <- stepwiseqtl(cross, method="hk", additive.only = T, penalties=c(lod_thresh, 0, 0))
# 
# 
# plot_qtl <- function(cross_file, blues_file) {
#   # c2 <- subset(cross_file, ind = (cross_file$pheno$genotype %in% blues_file$genotype))
#   # c2$pheno <- merge(c2$pheno[1], blues_file, by='genotype', all=F)
#   c2 <- jittermap(c2)
#   c2 <- calc.genoprob(c2, step=1)
#   print(nind(c2))
#   scan <- scanone(c2, pheno.col=grep("Height", colnames(c2$pheno)), method="hk")
#   plot(scan)
#   return(scan)
# }
# cross_old <- read.cross(format="csv",file="G:/My Drive/Nico_PhD.lnk/data/genotype/final_Gmap/UX1992_GMap.csv",
#                         estimate.map=FALSE, na.strings=c("-","NA"),
#                         genotypes=c("A", "H", "B"), crosstype = "riself")
# cross <- read.cross(format="csv",file=glue("linkage_map/maps/UX1992_linkage_map.csv"),
#                     estimate.map=FALSE, na.strings=c("-","NA"),
#                     genotypes=c("AA","H","BB"), crosstype="riself")
# cross_new <- read.cross(format="csv",file=glue("linkage_map/maps/UX1992.new_linkage_map.csv"),
#                         estimate.map=FALSE, na.strings=c("-","NA"),
#                         genotypes=c("A","H","B"), crosstype="riself")
# blues <- read.delim("data/blues.csv", sep=",") %>%
#   rename(genotype = Entry) %>%
#   select(-Cross_ID)
# oldblue <- read.delim("G:/My Drive/Nico_PhD.lnk/data/phenotype/2023_phenotype/trait_BLUEs.csv", sep=",") %>%
#   rename(genotype = Entry)
# 
# co2 <- plot_qtl(cross_old, oldblue)
# plot_qtl(cross, oldblue)
# plot_qtl(cross_new, oldblue)
# co <- plot_qtl(cross_old, blues)
# plot_qtl(cross, blues)
# cn <- plot_qtl(cross_new, blues)
# plot_qtl(final_map, blues)
# 
# ##check blues agreement
# blues_check <- merge(blues[,c("genotype", "Height")], oldblue[,c("genotype", "Height")], by="genotype")
# cor(blues_check$Height.x, blues_check$Height.y)
# plot(blues_check$Height.x, blues_check$Height.y)
# 
# ##check genotype agreement
# length(intersect(final_map$pheno$genotype, cross_old$pheno$genotype))
# # new_6A <- cross_new$geno[['6A']][['data']]
# new_6A <- final_map$geno[['6A']][['data']]
# old_6A <- cross_old$geno[['6A.2']][['data']]
# length(intersect(colnames(new_6A), colnames(old_6A)))
# 
# plot_new <- data.frame(name = colnames(new_6A), x = as.numeric(sapply(strsplit(colnames(new_6A), "_"), "[", 2)), y = 1)
# plot_old <- data.frame(name = colnames(old_6A), x = as.numeric(sapply(strsplit(colnames(old_6A), "_"), "[", 2)), y = 2)
# # max(c(plot_new$new, plot_old$old))
# plot_df <- rbind(plot_new, plot_old)
# plot(plot_df$x, plot_df$y)
# text(c(148e6, 148e6, 416e6, 416e6), c(1.2,1.8,1.2,1.8), labels=c("Rht-25", "Rht-25", "Rht-24", "Rht-24"), col = c('red', 'red', 'green', 'green'))
# text(63193160, 1.9, labels = "old peak", col = "blue")
# points(63193160, 1.99, col="blue", pch=19, cex=1.5)
# points(63193160, 1.01, col="blue", pch=19, cex=1.5)
# 
# np <- filter(plot_new, x > 6e7 & x < 7e7)
# op <- filter(plot_old, x > 6e7 & x < 7e7)
# nd <- new_6A[, np$name]
# od <- old_6A[, op$name]
# row.names(nd) <- final_map$pheno$genotype
# row.names(od) <- cross_old$pheno$genotype
# 
# nd <- nd[row.names(nd) %in% row.names(od),]
# od <- od[row.names(od) %in% row.names(nd),]
# 
# nd <- nd[order(rownames(nd)), ]
# od <- od[order(rownames(od)), ]
# 
# 
# 
# 
# reshape_for_plot <- function(df, label) {
#   df %>%
#     mutate(individual = rownames(df)) %>%
#     pivot_longer(-individual, names_to = "marker", values_to = "genotype") %>%
#     mutate(set = label)
# }
# 
# # Prepare both matrices
# plot_df1 <- reshape_for_plot(as.data.frame(od), "Old Map")
# plot_df2 <- reshape_for_plot(as.data.frame(nd), "New Map")
# 
# # Combine
# plot_data <- rbind(plot_df1, plot_df2)
# 
# # Set factor levels to maintain consistent ordering
# # plot_data$individual <- factor(plot_data$individual, levels = rev(rownames(nd)))
# plot_data$individual <- factor(plot_data$individual, levels = rev(rownames(nd[order(nd[,'S6A_69179565']),])))
# # plot_data$marker <- factor(plot_data$marker, levels = colnames(nd))  # assumes same marker order
# 
# # Plot with ggplot2
# p <- ggplot(plot_data, aes(x = marker, y = individual, color = factor(genotype))) +
#   geom_point(shape = 15, size = 2, na.rm = TRUE) +
#   scale_color_manual(values = c(`1` = "red", `2` = "green"), na.translate = FALSE) +
#   theme_minimal(base_size = 10) +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     panel.grid = element_blank(),
#     strip.text = element_text(size = 12)
#   ) +
#   facet_wrap(~set, scales = "free_x") +
#   labs(title = "Genotype Calls by Marker", x = "Marker", y = "Individual", color = "Genotype")
# 
# print(p)
