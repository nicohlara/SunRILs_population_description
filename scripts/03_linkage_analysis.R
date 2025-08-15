
library(here)
library(qtl)
library(qtl2)
library(dplyr)
library(glue)
library(data.table)
library(purrr)
library(ggplot2)
library(gaston)
library(stringr)
library(tidyr)
library(gridExtra)
library(scales)
source("scripts/linkage_map_helper_functions.R")

setwd(here())

pedigree <- read.csv("data/cross_info.csv", header=F, col.names = c("Cross_ID", "Parent_1", "Parent_2"))
chrom_lengths <- read.delim("data/chromosome_lengths.tsv")

chrom_lengths_split <-  read.delim("data/chrom_lengths_split.txt")



###note, UX2023 doesn't seem to have an alt parent (AGS2000)
#standardize input if necessary
for (fam in pedigree$Cross_ID) {
  print(fam)
  out_file <- glue("linkage_map/maps/{fam}_linkage_map")
  lines <- readLines(glue("{out_file}.csv"))
  genodata <- lines[4:length(lines)]
  genodata <- gsub("AA", "A", genodata)
  genodata <- gsub("BB", "B", genodata)
  ##make NAM founder parent always "A"
  geno_table <- read.csv(textConnection(genodata), header=F, stringsAsFactors = F)
  gt2 <- geno_table[1:6]
  parent_priority <- c("HILLIARD", "GA06493-13LE6", "NC08-23383")
  parent_present <- parent_priority[parent_priority %in% geno_table$V1][1]
  switched <- mapply(
    function(col_data, parent_allele) {
      if (parent_allele == "B") chartr("AB", "BA", col_data) else col_data
    },
    geno_table[, 7:ncol(geno_table)],
    as.character(geno_table[geno_table$V1 == parent_present, 7:ncol(geno_table), drop = FALSE])
  )
  ##save output
  genodata <- apply(cbind(gt2, switched), 1, paste, collapse = ",")
  
  writeLines(c(lines[1:3], genodata), con = glue("{out_file}_conv.csv"))
  print(length(lines))
  print(length(genodata)+3)
  ##remove monomorphic markers
  gt <- read.delim(glue("{out_file}_conv.csv"), sep=",")
  mt <- gt[grep("UX", gt$genotype, invert=T),]
  mt <- mt[!is.na(mt$index),]
  missing_one_parent <- which((mt[1, ] == "-" | mt[2, ] == "-") & !(mt[1, ] == "-" & mt[2, ] == "-"))
  mt_imputed <- mt
  for (col in missing_one_parent) {
    parent1 <- mt[1, col]
    parent2 <- mt[2, col]
    pop_calls <- gt[[col]][-c(1:2)]
    pop_alleles <- unique(pop_calls[pop_calls != "-"])
    if (length(pop_alleles) == 2) {
      if (parent1 == "-" && parent2 %in% c("A", "B")) {
        mt_imputed[1, col] <- ifelse(parent2 == "A", "B", "A")
      } else if (parent2 == "-" && parent1 %in% c("A", "B")) {
        mt_imputed[2, col] <- ifelse(parent1 == "A", "B", "A")
      }
    } else {
      mt_imputed[1, col] <- NA
      mt_imputed[2, col] <- NA
    }
  }
  keep_markers <- !is.na(mt_imputed[1, ]) & !is.na(mt_imputed[2, ])
  gt2 <- rbind(gt[c(1:2, grep("UX", gt$genotype)),keep_markers], mt_imputed)
  # mt2 <- gt[,!(mt[1,] == mt[2,] & mt[1,] != "-")]
  write.table(gt2, file = glue("{out_file}_conv.csv"), quote = F, row.names = F, col.names = T, sep=",", na="")
}


##note: UX2023 needs remaking
for (fam in pedigree$Cross_ID) {
  print(fam)
  cross <- read.cross(format="csv",file=glue("linkage_map/maps/{fam}_linkage_map_conv.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
  SunCross <- convert2cross2(cross)
  ##calculate QTL probabilities using geno prob
  pr <- calc_genoprob(SunCross, map = SunCross$gmap, cores=4)
  # bonf_threshold <- -log10((0.1 / totmar(cross)))
  ##calculate a kinship matrix of relationship among individuals
  ##Perform genome scan using Haley-Knott regression
  for (trait in colnames(SunCross$pheno)[-1]) {
    sig_threshold <- summary(scan1perm(pr, SunCross$pheno[, trait], n_perm=1000, cores=0), alpha=c(0.1))
    out <- scan1(pr, SunCross$pheno[, trait])
    pks <- find_peaks(out, SunCross$gmap, peakdrop=0.5, expand2markers=T, drop=0.5, threshold = sig_threshold)
    print(nrow(pks))
    if (nrow(pks) > 0 ) {
      Pos <- c()
      for (i in 1:nrow(pks)) {
        Pos <- c(Pos, find_marker(SunCross$gmap, pks[i,'chr'], pos=pks[i,'pos']))
      }
      pks$Pos <- Pos
      pks$cross <- fam; pks$trait <- trait
      if (exists('peak_df')) {peak_df <- rbind(peak_df, pks)} else {peak_df <- pks}
      for (chr in unique(pks$chr)) {
        effects <- data.frame(scan1coef(pr[,chr], SunCross$pheno[, trait] ))
        effects$Pos <- row.names(effects)
        effects$cross <- fam; effects$trait <- trait
        if ('NC08-23383' %in% SunCross$covar$genotype) {
          effects$NC08 <- cross$geno[[chr]]$data[grep("NC08", cross$pheno$genotype),]} else {effects$NC08 <- NA
          }
        if ('HILLIARD' %in% SunCross$covar$genotype) {
          effects$HILLIARD <- cross$geno[[chr]]$data[grep("HILLIARD", cross$pheno$genotype),]} else {effects$HILLIARD <- NA
          }
        if ('GA06493-13LE6' %in% SunCross$covar$genotype) {
          effects$GA13LE6 <- cross$geno[[chr]]$data[grep("GA06493-13LE6", cross$pheno$genotype),]} else {effects$GA13LE6 <- NA
          }
        if (exists('peak_effects')) {peak_effects <- rbind(peak_effects, effects)} else {peak_effects <- effects}
        
      }
    }
    png(filename=glue("figures/CIM_plots/{fam}_{trait}_CIM.png"),
        width=750*3, height=300*3, res=72*3,
        bg='transparent')
    par(mar=c(4,4,6,1))
    plot(out, SunCross$gmap, lodcolumn = 'pheno1', main=paste0("CIM of ", fam), bg='transparent')
    abline(h=sig_threshold, col='#CC0000')
    dev.off()
  }
}
peak_df <- merge(peak_df, peak_effects, by=c('Pos', 'cross', 'trait'))
write.table(peak_df, file="outputs/qtl2_cim.tsv", quote=F, sep="\t", row.names=F)
# peak_effects_subset <- peak_effects[paste0(peak_effects$Pos, peak_effects$cross, peak_effects$trait) %in% paste0(peaks$Pos, peaks$cross, peaks$trait),]
# write.table(peak_effects, file="outputs/qtl2_cim_peak_effects.tsv", quote=F, sep="\t", row.names=F)

## PLOT LOD, CHROMOSOME, AND LD HEATMAP
trait <- 'PM'
chr <- "1A"
chr_lod <- out[grep(chr, rownames(out)), ,drop=FALSE]
chr_map <- SunCross$gmap[[chr]]
lod_df <- tibble(
  marker = names(chr_map),
  cM = unname(chr_map),
  pos = as.numeric(sapply(strsplit(names(SunCross$gmap$'1A'), "_"), "[", 2)),
  lod = as.numeric(chr_lod)
)

gmat <- SunCross$geno[[chr]]
calc_ld_matrix <- function(mat) {
  mat <- scale(mat, center = TRUE, scale = FALSE)  # mean-center genotypes
  cor_mat <- cor(mat, use = "pairwise.complete.obs")
  ld_mat <- cor_mat^2
  return(ld_mat)
}
ld_matrix <- calc_ld_matrix(gmat)
markers <- lod_df$pos[order(lod_df$pos)]
ld_df <- expand.grid(
  marker1 = markers,
  marker2 = markers
) %>%
  filter(marker1 < marker2) %>%  # only unique pairs (lower triangle)
  mutate(
    index1 = match(marker1, markers),
    index2 = match(marker2, markers),
    offset = index2 - index1  # number of markers between
  ) %>%
  mutate(
    x = (marker1 + marker2) / 2,
    y = -offset,
    ld = mapply(
      function(m1, m2) ld_matrix[glue("S1A_{m1}"), glue("S1A_{m2}")],
      marker1,
      marker2))

p1 <- ggplot(lod_df, aes(x = pos, y = lod)) +
  geom_line(linewidth = 1.2, color = "#2c7fb8") +
  labs(y = "LOD", x = NULL, title = glue("QTL Scan of {trait} on {chr} in {fam}")) +
  theme_minimal(base_size = 12) +
  theme(plot.margin = margin(10, 10, 0, 10))
p2 <- ggplot(lod_df, aes(x = pos, y = 0)) +
  geom_segment(aes(xend = pos, yend = 0.3), color = "black") +
  geom_segment(aes(x = min(pos), xend = max(pos), y = 0, yend = 0), size = 1) +
  theme_void() +
  theme(plot.margin = margin(0, 10, 0, 10))
p3 <- ggplot(ld_df, aes(x = x, y = y, fill = ld)) +
  geom_point(shape = 21, size = 4, color =  "transparent") +
  scale_fill_gradient(low = "white", high = "red", name = "ld") +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )
p4 <- grid.arrange(p1, p2, p3, heights = c(3, 0.5, 5), ncol = 1)
ggsave(glue("figures/QTL_scan_{chr}_{trait}_{fam}.png"), p4)



##full heatmaps
for (fam in pedigree$Cross_ID) {
  print(fam)
  cross <- read.cross(format="csv",file=glue("linkage_map/maps/{fam}_linkage_map_conv.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
  png(filename=glue("figures/linkage_heatmaps/{fam}.png"), width=750, height=700)
  heatMap(cross, lmax=10)
  dev.off()
}


# png(filename="figures/UX1995_1A_PM.png",
#     width=750*3, height=300*3, res=72*3,
#     bg='transparent')
# par(mar=c(4,4,6,1))
# plot(out, chr='1A', SunCross$gmap, lodcolumn = 'pheno1', main="CIM of PM for UX1995 chrom 1A")
# abline(h=sig_threshold, col='#CC0000')
# dev.off()


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
  text(0.5,0.5, chr, cex=text_cex*1.5)
}

for (row in 1:nrow(ped)) {
  fam <- ped$Cross_ID[row]
  print(fam)
  plot.new()
  text(0,0.75, ped[row, 'Parent_2'], cex=text_cex*0.75, pos=4)
  text(0,0.25, ped[row, 'Parent_1'], cex=text_cex*0.75, pos=4)
  text(1, 0.5, ped[row, 'Cross_ID'], cex=text_cex, pos=2)
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


##get total markers/chrom for all linkage groups
stats <- data.frame(chrom = chrom_lengths$Chromosome)
stats_cM <- data.frame(chrom = chrom_lengths$Chromosome)

for (fam in ped$Cross_ID) {
  famstats <- read.delim(glue("linkage_map/maps/{fam}_stats.csv"), sep=" ")
  sums <- sapply(chrom_lengths$Chromosome, function(chrom) {
    subset_rows <- famstats[grep(chrom, rownames(famstats)),]
    sum(subset_rows$n.mar)
  })
  cM <-  sapply(chrom_lengths$Chromosome, function(chrom) {
    subset_rows <- famstats[grep(chrom, rownames(famstats)),]
    round(sum(subset_rows$length), 1)
  })
  stats[[fam]] <- sums
  stats_cM[[fam]] <- cM
  # famstats <- cbind( fam = rep(fam, nrow(famstats)), chrom = rownames(famstats), famstats)
}
stats <- rbind(stats, c('overall', colSums(stats[-1])))
stats_cM <- rbind(stats_cM, c('overall', colSums(stats_cM[-1])))

write.table(stats, file="outputs/linkage_map_marker_numbers.tsv", quote=F, sep="\t", row.names=F,)
write.table(stats_cM, file="outputs/linkage_map_cM_lengths.tsv", quote=F, sep="\t", row.names=F,)


## get bp positions of motifs
source("scripts/linkage_map_helper_functions.R")

chromosome <- '5B'
# fam <- 'UX2013'

crosses <- pedigree$Cross_ID
# crosses <- filter(pedigree, Parent_2 == "GA06493-13LE6")$Cross_ID

##plot chrom curves
for (fam in pedigree$Cross_ID) {
  cross <- read.cross(format="csv",file=glue("linkage_map/maps/{fam}_linkage_map.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("AA","H","BB"), crosstype="riself")
  for (chrom in head(rownames(summary.map(cross)), -1)) {
    plot_frame <- data.frame(cM = cross$geno[[chrom]]$map,
                             pos = sapply(strsplit(names(cross$geno[[chrom]]$map), "_"), "[", 2))
    png(filename = glue("figures/chrom_curves/{fam}_{chrom}_curve.png"),  
        width=400, height=400)
    plot(plot_frame$cM, plot_frame$pos, cex=2, pch = 20, xlab="cM", ylab="bp", main=glue("{fam} {chrom}"))
    dev.off()
  }
}



start_end <- c()

for (fam in crosses) {
  cross <- read.cross(format="csv",file=glue("linkage_map/maps/{fam}_linkage_map.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
  # plot_chrom_curve(cross, chromosome)
  # start_end <- c()
  for (chrom in grep(chromosome, rownames(summary.map(cross)), value=T)) {
    start_end <- c(start_end, lasso_marker_selector(cross, chrom))
  }
  # start_end
}
pos = as.numeric(sapply(strsplit(start_end, "_"), "[", 2))/1e6
min(pos); max(pos)


plot_frame <- data.frame(cM = cross$geno[[chromosome]]$map,
                         pos = sapply(strsplit(names(cross$geno[[chromosome]]$map), "_"), "[", 2))
par(cex.axis=1, mar=c(3,3,1,0))
png(filename = glue("figures/chrom_curves/{fam}_{chromosome}_curve.png"),  
    width=400, height=400)
plot(plot_frame$cM, plot_frame$pos, ylim=c(1, 1e9), cex=1, pch = 20, xlab="cM", main=glue("{fam} {chromosome}"), ylab="bp")
dev.off()

## generate initial monotonic maps from crosses
for (fam in pedigree$Cross_ID) {
  linkagemap <- read.delim(glue("linkage_map/maps/{fam}_linkage_map_conv.csv"), sep=",")
  linkagemap <- linkagemap[1:2, 7:ncol(linkagemap)]
  lm <- data.frame(chrom = as.character(linkagemap[1,]),
                   marker = colnames(linkagemap),
                   cM = round(as.numeric(linkagemap[2,]), 2),
                   pos = as.numeric(sapply(strsplit(colnames(linkagemap), "_"), "[", 2)))
  
  chroms_to_convert <- unique(sapply(strsplit(grep("\\.", unique(lm$chrom), value=T), "\\."), "[", 1))
  for (chr in chroms_to_convert) {
    add <- max(lm[lm$chrom == glue("{chr}.1"), "cM"]) + 500
    lm <- lm %>% mutate(cM = ifelse(chrom==glue("{chr}.2"), cM+add, cM))
  }
  lm <- lm %>% mutate(chrom = sapply(strsplit(chrom, "\\."), "[", 1))
  ##enforce cM/pos order
  lm_mono <- lm
  # diff <- 1
  # while (diff > 0) {
  #   r1 <- nrow(lm_mono)
  #   lm_mono <- lm_mono %>%
  #     arrange(chrom, cM) %>%
  #     group_by(chrom) %>%
  #     mutate(diffp = c(NA, diff(pos)),
  #            diffc = c(NA, diff(cM))) %>%
  #     filter(is.na(diffp) | is.na(diffc) | diffp > 0 | diffc > 0) %>%
  #     select(-c(diffp, diffc)) %>%
  #     ungroup()
  #   diff <- r1-nrow(lm_mono)
  #   print(glue("{r1}, {diff}"))
  # }
  write.table(lm_mono, glue("linkage_map/monotonic/{fam}_GBS_monotonic.map"), quote=F, sep="\t", row.names=F, col.names=F)
}


##add missing chromosomes in to monotonic maps
map_files <- list.files("linkage_map/monotonic", pattern = "_GBS_monotonic.map$", full.names = TRUE)
all_maps <- map_df(map_files, ~ {
  read_tsv(.x, col_names = c("chrom", "marker", "cM", "pos")) %>%
    mutate(source = basename(.x))
})

monotonic_consensus <- data.frame()
chroms <- paste0(rep(1:7, each=3), rep(c("A", "B", "D"), 7))
for (chr in chroms) {
  chr_filter <- filter(all_maps, chrom == chr)
  path <- marker_path_selector(select(chr_filter, cM, pos), chr)
  path$chrom <- chr
  monotonic_consensus <- rbind(monotonic_consensus, path)
}
monotonic_consensus <- mutate(monotonic_consensus, pos = round(pos,0)) %>%
  mutate(marker = glue("{chrom}_{pos}")) %>%
  select(chrom, marker, cM, pos)
write.table(monotonic_consensus, "linkage_map/monotonic/consensus_GBS_monotonic.map", quote=F, sep="\t", row.names=F, col.names=F)


##add missing chromosomes, filter for out of order markers
for (fam in pedigree$Cross_ID) {
  map <- read.table(glue("linkage_map/monotonic/{fam}_GBS_monotonic.map"), header=F, col.names= c("chrom", "marker", "cM", "pos"))
  for (chr in chroms) {
    if (!(chr %in% map$V1)) {
      add <- monotonic_consensus[monotonic_consensus$chrom == chr,]
      map <- rbind(map, add)
    }
  }
  
  lm_mono <- map 
  diff <- 1
  while (diff > 0) {
    r1 <- nrow(lm_mono)
    lm_mono <- lm_mono %>%
      arrange(chrom, cM) %>%
      group_by(chrom) %>%
      mutate(diffp = c(NA, diff(pos)),
             diffc = c(NA, diff(cM))) %>%
      filter(is.na(diffp) | is.na(diffc) | (diffp > 0 & diffc > 0)) %>%
      select(-c(diffp, diffc)) %>%
      ungroup()
    diff <- r1-nrow(lm_mono)
    print(glue("{r1}, {diff}"))
  }
  map <- lm_mono
  
  map <- map %>% mutate(pos = round(pos, 0)) %>%
    mutate(marker=glue("S{chrom}_{pos}"))
  
  
  write.table(map, glue("linkage_map/monotonic/{fam}_GBS_monotonic.map"), quote=F, sep="\t", row.names=F, col.names=F)
}
