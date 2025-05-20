
library(here)
library(qtl)
library(qtl2)
library(dplyr)

setwd(here())

pedigree <- read.csv("data/cross_info.csv")
chrom_lengths <- read.delim("data/chromosome_lengths.tsv")
# 
##standardize input if necessary
for (fam in pedigree$Cross_ID) {
  print(fam)
  out_file <- glue("linkage_map/maps/{fam}_linkage_map")
  lines <- readLines(glue("{out_file}.csv"))
  genodata <- lines[4:length(lines)]
  genodata <- gsub("AA", "A", genodata)
  genodata <- gsub("BB", "B", genodata)
  writeLines(c(lines[1:3], genodata), con = glue("{out_file}_conv.csv"))
  print(length(lines))
  print(length(genodata)+3)
}

for (fam in pedigree$Cross_ID) {
  print(fam)
  cross <- read.cross(format="csv",file=glue("linkage_map/maps/{fam}_linkage_map_conv.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
  SunCross <- convert2cross2(cross)
  ##calculate QTL probabilities using geno prob
  pr <- calc_genoprob(SunCross, map = SunCross$gmap, cores=4)
  bonf_threshold <- -log10((0.1 / totmar(cross)))
  ##calculate a kinship matrix of relationship among individuals
  ##Perform genome scan using Haley-Knott regression
  for (trait in colnames(SunCross$pheno)[2]) {
    out <- scan1(pr, SunCross$pheno[, trait])
    pks <- find_peaks(out, SunCross$gmap, peakdrop=0.5, expand2markers=T, drop=0.5, threshold = 1.5)
    if (nrow(pks) > 0 ) {
      Pos <- c()
      for (i in 1:nrow(pks)) {
        Pos <- c(Pos, find_marker(SunCross$gmap, pks[i,'chr'], pos=pks[i,'pos']))
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
    plot(out, SunCross$gmap, lodcolumn = 'pheno1', main=paste0("CIM of ", fam), bg='transparent')
    abline(h=bonf_threshold, col='#CC0000')
    dev.off()
  }
}
write.table(peaks, file="outputs/qtl2_cim.tsv", quote=F, sep="\t", row.names=F)

##create unified map of all linkage groups
chrom_df <- function(cross_map, chrom) {
  plot_frame <- data.frame(cM = cross_map$geno[[chrom]]$map,
                           pos = sapply(strsplit(names(cross_map$geno[[chrom]]$map), "_"), "[", 2))
  return(plot_frame)
}


plot_width <- nrow(chrom_lengths)+1
plot_depth <- nrow(pedigree)+1
plot_size <- 5
text_cex <- 2
plot_cex <- 1.5
res_val <- 25

ped <- pedigree %>% arrange(Parent_2, Parent_1)

png(filename = glue("figures/linkage_maps.png"),  
    width=plot_width*plot_size*res_val, height=plot_depth*plot_size*res_val)
par(mar=c(2,0,1,0))
plots <- plot_width*plot_depth
nf <- layout(matrix(c(1:plots), plot_depth, plot_width, byrow=T),
             heights=matrix(c(rep(plot_size, plot_width*plot_depth)), plot_depth, byrow=T),
             widths=matrix(c(rep(plot_size, plot_width*plot_depth)), plot_depth, byrow=T))
plot.new()

for (chr in chrom_lengths$Chromosome) {
  plot.new()
  text(0.5,0.75, chr, cex=text_cex)
}


for (row in 1:nrow(ped)) {
  fam <- ped$Cross_ID[row]
  print(fam)
  plot.new()
  text(0.5,0.75, ped[row, 'Parent_2'], cex=text_cex)
  text(0.5,0.25, ped[row, 'Parent_1'], cex=text_cex)
  cross <- read.cross(format="csv",file=glue("linkage_map/maps/{fam}_linkage_map.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
  groups <- summary.map(cross)
  for (chrom in chrom_lengths$Chromosome) {
    ymax <- chrom_lengths[chrom_lengths$Chromosome == chrom, 'End']*1.3
    grps <- grep(chrom, rownames(groups), value=T)
    if (length(grps) > 0) {  
      dfs <- lapply(grps, function(g) {
        d <- chrom_df(cross, g)
        d$grp <- g  # Keep track of group for later
        return(d)
      })
      if (length(dfs) > 1) {for (i in 2:length(dfs)) {m <- max(dfs[[i-1]][['cM']]); dfs[[i]][['cM']] <- dfs[[i]][['cM']] + m}}
      df_all <- do.call(rbind, dfs) %>%
        dplyr::mutate(cM = as.numeric(cM), pos = as.numeric(pos))
      if (max(df_all$pos) > ymax) {print(glue("too big: {chrom}"))}
      plot(df_all$cM, df_all$pos, cex=plot_cex, pch=20, xlab="cM", ylim=c(1, ymax), yaxt='n', ann=F)
      if (length(dfs) > 1) {for (i in 2:length(dfs)) {abline(v = min(dfs[[i]][['cM']]), col='red', cex=2)}}
    } else {
      plot.new()
    }
  }
}
dev.off()
