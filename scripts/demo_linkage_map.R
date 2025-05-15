### Demonstration linkage map creation for a cross object
## generally speaking, choosing two linkage groups per chromosome is not advised


##install needed packages
packages <- c("dplyr", "qtl", "ASMap", "stringr", "glue", "shiny", "plotly")
install.packages(setdiff(packages, rownames(installed.packages())))  

##load helper functions. Replace with your path to this file
source("scripts/linkage_map_helper_functions.R") 

##load genetic data
fam <- "UX1989"
## if using gaston to read in vcf, convert to matrix and then to cross
# vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
# geno <- as.data.frame(as.matrix(gaston_object), stringAsFactors = F) 
# cross <- convert_vcf_to_cross(geno, blues)
# write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")

###filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)

##change based on histogram
missing_threshold <- 32500
abline(v=missing_threshold, col='red')
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < missing_threshold)

###filter by marker stats
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
##Aiming for ~4k markers with >100 markers per chromosome
##change value based on plot, will take any marker with LESS than miss_threshold missing data.
miss_threshold = 0.65
abline(v=miss_threshold, col='red')
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross1$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25, display.markers=F)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
###check on marker numbers before proceeding

##Create linkage maps. Click on any linkage groups that look good
##TIP: PULL OPENED TAB TO ITS OWN WINDOW! 
##    Marker selection opens tons of new tabs and its nice to delete them all at once
##optional min_marker parameter: default is 10, any linkage groups with less are deleted
map1 <- make_map(SunCross4, p_value = 1e-12)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)

##optional: create a second map with the remaining chromosome groups at a higher p-value
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- make_map(map2, p_value = 1e-15)
keep2 <- linkage_group_selector(map2)
final_map2 <- subset(map2, keep2)
final_map <- combineMap(final_map, final_map2, id='genotype')

## finalize map, rename linkage groups, join or flip chroms, rogue out markers etc.
## flips chromosomes
to_flip <- linkage_group_selector(final_map)
final_map2 <- flip.order(final_map, c(to_flip))

## join linkage groups.
##NEED TO SELECT IN ORDER WITHIN EACH CHROMOSOMAL GROUP
## e.g. if the plots are 2D long arm, 2D short arm, select 2Ds THEN 2Dl or they will be stitched backwards
## often results in VERY long chromosomes
filter <- TRUE
max_distance <- 500
to_join <- linkage_group_selector(final_map)
final_map3 <- final_map2
for (grp in unique(gsub("\\..*", "",to_join))) {
  L <- list()
  L[[grp]] <-  grep(grp, to_join, value=T)
  print(L)
  final_map4 <- mergeCross(final_map3, merge = L)
  if (filter == TRUE) {
    distance <- summary.map(quickEst(final_map4, chr = grp, map.function = "kosambi"))[grp, 'length']
    if (distance < max_distance) {final_map3 <- final_map4}
  } else {
    final_map3 <- final_map4
  }
}
summary.map(final_map3)
##if this looks good, proceed

##clean out any rogue markers in linkage groups
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

##select final linkage groups
final_map5 <- linkage_group_selector(final_map4)
final_map <- subset(final_map4, final_map5)

##re-estimate distances, print stats, and save plots, stats, and the map
output_directory <- "linkage_map/maps/"
## this gives pretty bad results with stitched linkage groups
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map))
write.table(summary.map(final_map), file = glue("{output_directory}{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path=glue("{output_directory}"))
write.cross(final_map, "csv", filestem=glue("{output_directory}{fam}_linkage_map"))
