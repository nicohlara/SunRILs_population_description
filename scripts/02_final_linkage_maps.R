# library("dplyr")
# library("qtl")
# library("ASMap")
library("gaston")
# library(stringr)
# library(glue)
# library(shiny)
# library(plotly)
library(here)

setwd(here())

source("scripts/linkage_map_helper_functions.R") 
###notes:
##selection apps should open in browser. When done press 'Done' at bottom of page, then page can be exited and results used
##Phenotype needs to be incorporated BEFORE creating a cross object or it doesn't seem to line up properly

### Coerce all vcf files to cross objects
pedigree <- read.delim("data/cross_info.csv", sep=",")
blues <- read.delim("data/blues.csv", sep=",") %>%
  rename(genotype = Entry) %>%
  select(-Cross_ID)

redone_fams <- c("UX1992")
# for (fam in pedigree$Cross_ID) { ##original
for (fam in redone_fams) {
  print(fam)
  # vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE)
  vcf <- read.vcf(glue("data/processed_vcf/SunRILs_raw_{fam}_subset_filt.vcf.gz"), convert.chr = F)
  if (fam == "UX2031") {
    vcf@ped$id <- gsub("UX2031-99-", "UX2031-99", unique(vcf@ped$id))
  }
  genotype <- as.data.frame(as.matrix(vcf), stringAsFactors = F) 
  cross <- convert_vcf_to_cross(genotype, blues)
  # cross <- convert_vcf_to_cross(vcf, blues)
  write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)
}

filter <- data.frame(popmiss = c(30, 30, 22.5, 27.5, 20, 25, 30, 25, 25, 27.5, 27.5, 30, 25, 25, 35),
                     markmiss = c(0.5, .5, .65, .6, .55, .65, .4, .5, .55, .8, .75, .7, .7, .75, .7))
rownames(filter) <- pedigree$Cross_ID

redo_filter <- data.frame(popmiss = c(35),
                          markmiss= c(.75))
rownames(redo_filter) <- redone_fams
##filter down the markers
# for (fam in pedigree$Cross_ID) {
for (fam in redone_fams) {
  # pop_missing <- filter[fam, 'popmiss']*1000
  # miss_threshold <- filter[fam, 'markmiss']
  pop_missing <- redo_filter[fam, 'popmiss']*1000
  miss_threshold <- redo_filter[fam, 'markmiss']
  SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                        estimate.map=FALSE, na.strings=c("-","NA"),
                        genotypes=c("A","H","B"), crosstype="riself")
  ##filter by individuals
  pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
  # hist(pg$stat$miss)
  SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < pop_missing)
  ##filter by markers
  nt.bymar <- ntyped(SunCross1, 'mar')
  # hist(nt.bymar/length(SunCross1$pheno$genotype))
  ###Aiming for ~4k markers with >100 markers per chromosome
  SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
  SunCross3<-pullCross(SunCross2,type="co.located")
  pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                      "bonf", type = "l", cex = 0.25)
  SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
  
  print(fam)
  print(c(totmar(SunCross), totmar(SunCross1), totmar(SunCross2), totmar(SunCross3), totmar(SunCross4)))
  print("")
  write.cross(SunCross4, "csv", filestem=glue("linkage_map/filt_premap/{fam}_filtered"))
}

# power <- data.frame(p1 = c(10, 14, 12, 10, 10, 10, 10, 10, 6, 10, 10, 10, 10, 10, 10),
                    # p2 = c(NA, 22, NA, NA, 15, 15, NA, NA, 15, 15, NA, NA, 25, NA, NA))
# rownames(power) <- pedigree$Cross_ID
stats <- list()

##make maps
# for (fam in pedigree$Cross_ID[-c(1:7)]) {
for (fam in redone_fams) {
  print(fam)
  SunCross<- read.cross(format="csv",file="linkage_map/SunRILs_filtered.csv",
                        estimate.map=FALSE, na.strings=c("-","NA"),
                        genotypes=c("AA","H","BB"), crosstype="riself")
  # SunCross <-  read.cross(format="csv",file=glue("linkage_map/filt_premap/{fam}_filtered.csv"),
  #                         estimate.map=FALSE, na.strings=c("-","NA"),
  #                         genotypes=c("AA","H","BB"), crosstype="riself")
  # ##create map
  print("Choose chromosomes")
  p1 <- as.numeric(readline("First p-value? 10^-n "))
  map1 <- make_map(SunCross, p_value = 10^-p1)
  keep <- linkage_group_selector(map1)
  final_map <- subset(map1, keep)

  summaryMap <- summary.map(map1)
  remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
  
  if (length(remove) > 1) {
    ##optional: create a second map with the remaining chromosome groups at a higher p-value
    p2 <-  as.numeric(readline(glue("Second p-value? 10^-n, first was {p1}")))
    print("Choose chromosomes pt 2")
    map2 <- subset(SunCross, chr=unique(head(gsub("\\..*", "", remove), -1)))
    map2 <- make_map(map2, p_value = 10^-p2)
    keep2 <- linkage_group_selector(map2)
    final_map2 <- subset(map2, keep2)
    final_map2$pheno <- final_map2$pheno[1]
    final_map <- combineMap(final_map, final_map2, id='genotype')
  }
  
  ## finalize map, rename linkage groups, join or flip chroms, rogue out markers etc.
  print("Choose which chrom to flip")
  to_flip <- linkage_group_selector(final_map, zoom=T)
  if (length(to_flip) > 0) {final_map2 <- flip.order(final_map, c(to_flip))} else {final_map2 <- final_map}
  
  ## join linkage groups.
  filter <- TRUE
  max_distance <- 500
  print("choose which linkage group to combine")
  to_join <- linkage_group_selector(final_map2)
  final_map3 <- final_map2
  for (grp in unique(gsub("\\..*", "",to_join))) {
    L <- list()
    L[[grp]] <-  grep(grp, to_join, value=T)
    print(L)
    final_map4 <- mergeCross(final_map3, merge = L)
    if (filter == TRUE) {
      distance <- summary.map(quickEst(final_map4, chr = grp, map.function = "kosambi"))[grp, 'length']
      print(distance)
      if (distance < max_distance) {final_map3 <- final_map4}
    } else {
      final_map3 <- final_map4
    }
  }
  summary.map(final_map3)
  ##clean out any rogue markers in linkage groups
  print("clean out markers")
  markers_to_remove <- c()
  for (chr in chrnames(final_map3)) {
    message("Selecting markers to remove for chromosome: ", chr)
    selected <- lasso_marker_selector(final_map3, chr)
    if (!is.null(selected)) {
      markers_to_remove <- c(markers_to_remove, selected)
    }
  }
  length(markers_to_remove)
  final_map4 <- drop.markers(final_map3, markers_to_remove)
  ##select final linkage groups in order
  print("finalize selection")
  final_chroms <- linkage_group_selector(final_map4)
  final_map5 <- subset(final_map4, final_chroms)
  final_map5$geno <- final_map5$geno[final_chroms]
  
  df <- data.frame(original = final_chroms) %>%
    mutate(base = str_extract(original, "^[^\\.]+")) %>%  # extract base before dot
    group_by(base) %>%
    mutate(new = if (n() == 1) base else paste0(base, ".", row_number())) %>%
    ungroup()
  names(final_map5$geno) <- df$new
  ##re-estimate distances, print stats, and save plots, stats, and the map
  output_directory <- "linkage_map/maps/"
  ## this gives pretty bad results with stitched linkage groups
  final_map <- quickEst(final_map5, map.function = "kosambi")
  print(summary.map(final_map))
  write.table(summary.map(final_map), file = glue("{output_directory}{fam}_stats.csv"), quote=F, row.names=T)
  plot_linkage_map(final_map, fam, file_path=glue("{output_directory}"))
  write.cross(final_map, "csv", filestem=glue("{output_directory}{fam}_linkage_map"))
  
  stats[[glue("{fam}_p1")]] <- p1
  stats[[glue("{fam}_keep1")]] <- c(keep)
  if (length(remove) > 1) {stats[[glue("{fam}_keep2")]] <- keep2; stats[[glue("{fam}_p2")]] <- p2}
  stats[[glue("{fam}_flipped")]] <- to_flip
  stats[[glue("{fam}_joined")]] <- to_join
  stats[[glue("{fam}_removed_markers")]] <- markers_to_remove
  stats[[glue("{fam}_final_groups")]] <- final_chroms
  
  
}

sink(glue("{output_directory}/linkage_creation_details.txt"))
print(stats)
sink()