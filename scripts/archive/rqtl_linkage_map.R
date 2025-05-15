library(tidyverse)
library("qtl")
library("ASMap")
library("gaston")
# library(bwardr)
library(glue)

convert_vcf_to_cross <- function(gaston_object) {
  geno <- as.data.frame(as.matrix(gaston_object), stringAsFactors = F)
  geno[geno==0] <- "A"
  geno[geno==1] <- "H"
  geno[geno==2] <- "B"
  geno[is.na(geno)] <- "-"
  geno <- cbind(1:dim(geno)[1], geno)
  colnames(geno)[1] <- "index"
  geno <- rbind(c("",gaston_object@snps$chr), geno)
  rownames(geno)[1] <- ""
  geno <- cbind(rownames(geno), geno)
  colnames(geno)[1] <- 'genotype'
  return(geno)
}  
plot_chrom_curve <- function(cross_map, chrom) {
  plot_frame <- data.frame(cM = cross_map$geno[[chrom]]$map,
                           pos = sapply(strsplit(names(cross_map$geno[[chrom]]$map), "_"), "[", 2))
  par(cex.axis=1.5)
  plot(plot_frame$cM, plot_frame$pos, ylim=c(1, 1e9))
}

plot_linkage_map <- function(cross_object, family, file_path="linkage_map/cross_objects/") {
  summaryMap <- summary.map(cross_object)
  plot_width <- nrow(summaryMap)-1
  plot_depth <- 2
  plot_size <- 4
  text_cex <- 4
  res_val <- 1
  pdf(file=glue("{file_path}/{family}_plots.pdf"), width=plot_width*plot_size*res_val, height=plot_depth*plot_size*res_val)
  nf <- layout(matrix(c(1:(plot_width*plot_depth)), nrow = plot_depth, ncol = plot_width, byrow=F), 
               heights=matrix(c(rep(plot_size/4,plot_width), rep(plot_size, plot_width*(plot_depth-1))), plot_depth, byrow=T), 
               widths=matrix(rep(c(rep(plot_size, plot_width)), plot_depth), nrow=plot_depth))
  par(mar=c(2,0,1,0))
  
  for (chrom in row.names(summaryMap)[-nrow(summaryMap)]) { 
    plot.new()
    text(0.5,0.75, chrom, cex=text_cex)
    text(0.5,0.25, summaryMap[chrom, 'n.mar'], cex=text_cex-1)
    plot_chrom_curve(cross_object, chrom)
  }
  dev.off()
}

remake_map <- function(map, p_value) {
  map2 <- mstmap(map, id='genotype', bychr = T,
                 anchor = T,
                 dist.fun='kosambi',
                 objective.fun = 'COUNT',
                 p.value = p_value,
                 noMap.dist = 15,
                 noMap.size = 0,
                 miss.thresh = 1,
                 mvest.bc=F,
                 detectBadData=F,
                 return.imputed=F,
                 as.cross=T, pop.type = 'RIL6',
                 trace=T)
  summaryMap <- summary.map(map2)
  ##remove small/incomplete chromosomes
  chrnamesDROP<-rownames(subset(summaryMap,summaryMap$n.mar < 10))
  markernamesDROP<-markernames(map2,chr=chrnamesDROP)
  map2<-drop.markers(map2,markers=markernamesDROP)
  summaryMap <- summary.map(map2)
  print(summaryMap)
  plot_linkage_map(map2, fam)
  return(map2)
}

### ------------------------------------------------------------------------------------ ###
fam <- "UX1989"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.30
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-8
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-4)
keep <- c('1A.2', '1B.4', '1D.2', '2A.2', '2B.2', '2D.4', '3A.3', '3B.2', '3D.2', '4A.2', '4B.2', '5B.2', '6A.1', '6B.2', '6D.4', '7B.1', '7D.6' )
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-8)
keep <- c('5A.4', '7A.6')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]


##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX2023"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.15
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-5
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-5)
keep <- c('1D.6', '2D.1', '3A.2', '3B.3', '3D.2', '4A.3', '4B.1', '5A.1', '5B.2', '5D.12', '6A.2', '6D.2', '7D.10')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-10)
keep <- c('1A.2', '1B.5', '6B.1', '7A.5', '7B.2')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map3 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-15)
keep <- c('2B.5')
final_map3 <- subset(map3, keep)
final_map <- combineMap(final_map, final_map3, id='genotype')
summaryMap <- summary.map(map3)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map4 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-16)
keep <- c('2A.9')
final_map4 <- subset(map4, keep)
final_map <- combineMap(final_map, final_map4, id='genotype')
# summaryMap <- summary.map(map4)
# remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]


##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX1991"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE)
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.30
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-8
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-4)
keep <- c('1A.1', '1D.2', '3A.1', '3B.2', '3D.4', '4B.4', '5A.2', '5D.4', '6A.1', '6B.2', '6D.2', '7B.2', '7D.1')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-7)
keep <- c('1B.2', '2A.1', '2B.1', '2D.3', '4A.2' )
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map3 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-15)
keep <- c('5B.3', '7A.5')
final_map3 <- subset(map3, keep)
final_map <- combineMap(final_map, final_map3, id='genotype')
# summaryMap <- summary.map(map3)
# remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c('5D'))
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX1992"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.10
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
# segregation_threshold = 1e-5
# SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
# SunCross2
SunCross3<-pullCross(SunCross1,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-2)
# keep <- c('1A.2', '1B.1', '2B.1', '2D.1', '4B.3', '5B.2', '6A.1', '6D.2', '7D.1')
keep <- c('1D.3', '2A.2', '2D.2', '3B.2', '4A.3', '5A.1', '5B.1', '6B.3', '6D.2', '7D.2')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-5)
# keep <- c('1D.2', '2A.1', '3A.2', '3B.1', '4A.2', '5A.2', '6B.1', '7A.2', '7B.2')
keep <- c("1A.3", '1B.4', '2B.3', '3A.4', '4B.5', '6A.5', '7A.1', '7B.2')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX1993"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.35
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-10
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-5)
keep <- c('1B.1', '1D.4', '2A.3', '2D.3', '4A.1', '4B.2', '5A.1', '5D.7', '6A.2', '6D.2', '7A.1', '7B.5' )
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-9)
keep <- c('1A.1', '3A.4', '3B.2', '5B.4', '6B.3', '7D.1')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map3 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-15)
keep <- c('2B.2')
final_map3 <- subset(map3, keep)
final_map <- combineMap(final_map, final_map3, id='genotype')
# summaryMap <- summary.map(map3)
# remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX1994"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.45
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-8
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-3)
keep <- c('1A.4', '1B.2', '1D.1', '2A.1', '2B.1', '2D.1', '3A.1', '3B.1', '3D.3', '4A.2', '4B.3', '4D.1', '5A.3', '5B.1', '5D.8', '6D.4', '7A.3', '7B.1')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- subset(SunCross3, chr=head(gsub("\\..*", "", remove), -1))
map2 <- remake_map(map2, 1e-6)
keep <- c('6A.1')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')

##check and polish map
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c('2D.1' ))
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX1995"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.25
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-5
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-2)
keep <- c('1A.2', '1D.1', '2A.1', '2B.2', '2D.3', '3D.1', '4D.1', '6A.1' )
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- remake_map(map2, 1e-5)
keep <- c('1B.1', '3A.3', '3B.3', '4A.2', '4B.2', '5A.1', '5B.1', '7B.2', '7D.8' )
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map3 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map3 <- remake_map(map3, 1e-9)
keep <- c('6B.2', '7A.1')
final_map3 <- subset(map3, keep)
final_map <- combineMap(final_map, final_map3, id='genotype')
# summaryMap <- summary.map(map3)
# remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX1997"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.40
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-9
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-4)
keep <- c('1A.2', '1D.1', '3D.3', '4A.1', '4B.1', '4D.2', '5D.1', '6D.2', '7D.4')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- remake_map(map2, 1e-8)
keep <- c('1B.2', '2B.2', '2D.1', '5A.1', '5B.3', '6A.2', '6B.2', '7A.1', '7B.2')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map3 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map3 <- remake_map(map3, 1e-12)
keep <- c('3A.5', '3B.2')
final_map3 <- subset(map3, keep)
final_map <- combineMap(final_map, final_map3, id='genotype')
# summaryMap <- summary.map(map3)
# remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX2000"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.25
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-5
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-2)
keep <- c('1B.1', '1D.2', '2D.2', '4B.2', '4D.2', '7D.3')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- remake_map(map2, 1e-5)
keep <- c('2B.1', '3A.2', '3B.3', '3D.1', '4A.1', '5A.2', '5B.1', '6A.2', '6D.2')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map3 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map3 <- remake_map(map3, 1e-8)
keep <- c('1A.1', '2A.2', '6B.1', '7A.2', '7B.2')
final_map3 <- subset(map3, keep)
final_map <- combineMap(final_map, final_map3, id='genotype')
# summaryMap <- summary.map(map3)
# remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c('3D'))
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX2010"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.25
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-9
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-6)
keep <- c('1A.1', '1B.5', '1D.3','2D.3', '3A.2', '6A.5',  '6D.5', '7D.4')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- remake_map(map2, 1e-12)
keep <- c('2B.7', '4A.7', '4B.3', '5A.5', '5B.4', '6B.8', '7A.9', '7B.2')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map3 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map3 <- remake_map(map3, 1e-16)
keep <- c('2A.2', '3B.12')
final_map3 <- subset(map3, keep)
final_map <- combineMap(final_map, final_map3, id='genotype')
# summaryMap <- summary.map(map3)
# remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX2012"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.20
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-5
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-4)
keep <- c('1A.1', '1B.3', '1D.4', '2A.2', '2D.1', '3A.2', '4A.2', '4B.2', '6A.1', '6D.2', '7A.4', '7B.1', '7D.4' )
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- remake_map(map2, 1e-8)
keep <- c('2B.2', '3B.2', '5A.1', '5B.2',  '6B.3')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX2013"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.25
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-5
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-3)
keep <- c('1B.2', '1D.1', '2A.1', '2B.2', '3B.2', '3D.1', '4B.2', '5B.1', '5D.2', '6A.2', '6B.1', '6D.3')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- remake_map(map2, 1e-6)
keep <- c('1A.1', '2D.3', '3A.3', '4A.1', '5A.2', '7A.1', '7B.3' )
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX2026"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.25
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-5
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-2)
keep <- c('1A.1', '1D.2', '2D.2', '3A.2', '4B.2', '5D.1')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-4)
keep <- c('1B.1', '3D.3', '5A.2', '6A.3', '6D.2', '7B.1', '7D.3')
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map3 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-12)
keep <- c('2A.1', '2B.1', '3B.3', '4A.2', '5B.1', '7A.1')
final_map3 <- subset(map3, keep)
final_map <- combineMap(final_map, final_map3, id='genotype')
# summaryMap <- summary.map(map3)
# remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX2029"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.25
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-5
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-3)
keep <- c('1A.2', '1B.1', '1D.3', '2A.1', '2D.1', '3A.1', '5A.1', '5B.1', '5D.5', '6A.2', '6B.1', '6D.3')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-6)
keep <- c('2B.2', '3B.1', '3D.4', '4A.3', '7A.1', '7B.2', '7D.6' )
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------------------ ###
fam <- "UX2031"
print(fam)
vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
cross <- convert_vcf_to_cross(vcf)
write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)

SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##remove markers with too much missing data
nt.bymar <- ntyped(SunCross, 'mar')
hist(nt.bymar/length(SunCross$pheno$genotype))
###edit this, aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.4
SunCross1 <- drop.markers(SunCross, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(1-miss_threshold)]))
SunCross1
segregation_threshold = 1e-5
SunCross2<-pullCross(SunCross1,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))
SunCross2
SunCross3<-pullCross(SunCross2,type="co.located")
SunCross3

##create linkage map
map1 <- remake_map(SunCross3, 1e-6)
keep <- c('1B.2', '1D.5', '2A.1', '3A.2', '3D.1', '4B.2', '4D.3', '5A.2', '5B.1', '6A.2', '6B.2', '6D.2', '7D.4')
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

map2 <- remake_map(subset(SunCross3, chr=unique(head(gsub("\\..*", "", remove), -1))), 1e-12)
keep <- c('1A.2', '2B.1', '2D.3', '3B.2', '4A.1', '5D.9', '7A.1', '7B.1' )
final_map2 <- subset(map2, keep)
final_map <- combineMap(final_map, final_map2, id='genotype')
summaryMap <- summary.map(map2)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]

##check and polish map
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
plot_linkage_map(final_map, fam)
final_map <- flip.order(final_map, c())
write.cross(final_map, "csv", filestem=glue("linkage_map/cross_objects/{fam}_linkage_map"))
### ------------------------------------------------------------------------------------ ###


################################################################################
## refine maps by removing markers
library(shiny)
library(plotly)

lasso_marker_selector <- function(cross, chr) {
  geno_pos <- pull.map(cross, chr=chr)[[1]]
  markers <- names(geno_pos)
  df <- data.frame(
    marker = markers,
    cM = geno_pos,
    pos = as.numeric(sapply(strsplit(markers, "_"), "[", 2))  # dummy y values just for plotting
  )
  # Create a temporary environment to hold result
  env <- new.env()
  env$selected_markers <- NULL
  
  ui <- fluidPage(
    h4(paste("Select markers to remove on chromosome", chr)),
    plotlyOutput("plot"),
    actionButton("done", "Done")
  )
  
  server <- function(input, output, session) {
    output$plot <- renderPlotly({
      plot_ly(df, x = ~cM, y = ~pos, type = "scatter", mode = "markers",
              text = ~marker, source = "select") %>%
        layout(dragmode = "lasso")
    })
    
    observeEvent(input$done, {
      selected <- event_data("plotly_selected", source = "select")
      if (!is.null(selected)) {
        env$selected_markers <- df$marker[selected$pointNumber + 1]
      }
      stopApp()
    })
  }
  
  runApp(shinyApp(ui, server), launch.browser = TRUE)
  
  return(env$selected_markers)
}

for (fam in cross_info$Cross_ID) {  
  print(fam)
  SunCross <- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_linkage_map.csv"),
                         estimate.map=FALSE, na.strings=c("-","NA"),
                         genotypes=c("A","H","B"), crosstype="riself")
  markers_to_remove <- c()
  for (chr in chrnames(SunCross)) {
    message("Selecting markers to remove for chromosome: ", chr)
    selected <- lasso_marker_selector(SunCross, chr)
    if (!is.null(selected)) {
      markers_to_remove <- c(markers_to_remove, selected)
    }
  }
  length(markers_to_remove)
  
  SunCross <- drop.markers(SunCross, markers_to_remove)
  SunCross <- quickEst(SunCross, map.function = "kosambi")
  print(summary.map(SunCross)['overall',])
  write.table(summary.map(SunCross), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
  plot_linkage_map(SunCross, fam, file_path="linkage_map/maps/")
  write.cross(SunCross, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
}
