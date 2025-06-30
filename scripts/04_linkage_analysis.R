
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
source("scripts/linkage_map_helper_functions.R")

setwd(here())

pedigree <- read.csv("data/cross_info.csv", header=F, col.names = c("Cross_ID", "Parent_1", "Parent_2"))
chrom_lengths <- read.delim("data/chromosome_lengths.tsv")

chrom_lengths_split <-  read.delim("data/chrom_lengths_split.txt")
# 
##standardize input if necessary
# for (fam in pedigree$Cross_ID) {
#   print(fam)
#   out_file <- glue("linkage_map/maps/{fam}_linkage_map")
#   lines <- readLines(glue("{out_file}.csv"))
#   genodata <- lines[4:length(lines)]
#   genodata <- gsub("AA", "A", genodata)
#   genodata <- gsub("BB", "B", genodata)
#   writeLines(c(lines[1:3], genodata), con = glue("{out_file}_conv.csv"))
#   print(length(lines))
#   print(length(genodata)+3)
# }

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
  for (trait in colnames(SunCross$pheno)[-1]) {
    out <- scan1(pr, SunCross$pheno[, trait])
    pks <- find_peaks(out, SunCross$gmap, peakdrop=0.5, expand2markers=T, drop=0.5, threshold = 1.5)
    print(nrow(pks))
    if (nrow(pks) > 0 ) {
      Pos <- c()
      for (i in 1:nrow(pks)) {
        Pos <- c(Pos, find_marker(SunCross$gmap, pks[i,'chr'], pos=pks[i,'pos']))
      }
      pks$Pos <- Pos
      pks$cross <- fam; pks$trait <- trait
      if (exists('peaks')) {peaks <- rbind(peaks, pks)} else {peaks <- pks}
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
peaks <- merge(peaks, peak_effects, by=c('Pos', 'cross', 'trait'))
write.table(peaks, file="outputs/qtl2_cim.tsv", quote=F, sep="\t", row.names=F)
# peak_effects_subset <- peak_effects[paste0(peak_effects$Pos, peak_effects$cross, peak_effects$trait) %in% paste0(peaks$Pos, peaks$cross, peaks$trait),]
# write.table(peak_effects, file="outputs/qtl2_cim_peak_effects.tsv", quote=F, sep="\t", row.names=F)




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


## heatmap of marker effects
pks <- read.delim("outputs/qtl2_cim.tsv") %>%
  filter(trait == 'PM') %>%
  mutate(phys_pos = round(as.numeric(sapply(strsplit(Pos, "_"), "[", 2))/1e6, 1),
         chr = gsub("\\..", "", chr)) %>%
  rename(chromosome = chr) %>%
  arrange(chromosome, phys_pos) %>%
  rename(Cross_ID = cross)

##divide into S/C/L
cls <- chrom_lengths_split
pks$chromosome <- apply(pks, 1, function(x) {y <- as.numeric(x[[17]])*1e6; cls[cls$Chromosome == x[[6]] & cls$Start <= y
                                                  & cls$End > y,'Arm']})

##add in GWAS summary
gwas_results <- read.delim("outputs/combined_GWAS.tsv") %>%
  mutate(lod = -log10(P.value))
gwas_results$chromosome <- apply(gwas_results, 1, function(x) {y <- as.numeric(x[[2]]); cls[cls$Chromosome == x[[1]] & cls$Start <= y
                                                                                            & cls$End > y,'Arm']})
gwas_results <- select(gwas_results, c(Position, Chromosome, chromosome, trait, lod)) %>%
  filter(trait == 'PM') %>%
  mutate(Cross_ID = 'GWAS', marker = glue("S{Chromosome}_{Position}") )
genotype <- read.vcf("data/SunRILs_prod_filt_imp.vcf.gz", convert.chr=F)
genotype@ped$family <- ifelse(grepl("UX", genotype@ped$id), str_sub(genotype@ped$id, 1, 6), 'Parent')
genotype <- select.snps(genotype, id %in% gwas_results$marker)
blues <- read.delim("data/blues.csv", sep=",") %>%
  mutate(Cross_ID = as.factor(Cross_ID)) %>% 
  rename(id = Entry)
genotype <- select.inds(genotype, id %in% blues$id)

geno <- cbind(id = genotype@ped$id, as.matrix(genotype)) %>%
  data.frame() 
geno_long <-data.frame(geno) %>%
  inner_join(blues[,c('id', 'PM')], by = "id") %>%
  pivot_longer(cols = -c(id, PM), names_to = "marker", values_to = "genotype")

# 2. Compute means per marker/genotype
geno_summary <- geno_long %>%
  group_by(marker, genotype) %>%
  summarise(mean_blue = mean(PM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = mean_blue, names_prefix = "mean_pheno_")

# View the result
geno_centered <- geno_summary %>%
  mutate(
    midpoint = (mean_pheno_0 + mean_pheno_2) / 2,
    adj_pheno_0 = mean_pheno_0 - midpoint,
    adj_pheno_2 = mean_pheno_2 - midpoint
  ) %>%
  select(marker, adj_pheno_0, adj_pheno_2)

gwas <- merge(gwas_results, geno_centered, by= 'marker', all=T) %>%
  rename(Pos = marker, AA = adj_pheno_0, BB = adj_pheno_2)

pks1 <- merge(pks, gwas, by=intersect(names(pks), names(gwas)), all=T)

peak_summary <- as.data.table(pks1)[, {
  idx <- which.max(abs(AA))
  .(
    effect_peak = pos[idx],
    phys_pos = phys_pos[idx],
    LOD = lod[idx],
    AA = AA[idx],
    BB = BB[idx],
    NC08 = NC08[idx],
    HILLIARD = HILLIARD[idx],
    GA13LE6 = GA13LE6[idx],
    nmar = uniqueN(pos)
  )
}, by = .(trait, chromosome, Cross_ID)]

peak_summary1 <- merge(peak_summary, pedigree, by='Cross_ID') %>%
  mutate(cross_label = paste(Parent_2, "-", Parent_1))
         # peak_label = paste0(chromosome, ":", peak_start, "-", peak_end)) 

##convert to NAM founder allele coding
peak_summary1 <- peak_summary1 %>%
  rename("NC08-23383" = NC08,
         "GA06493-13LE6" = GA13LE6) %>% 
  rowwise() %>%
  mutate(
    genotype_val = get(Parent_2),
    missing_geno = ifelse(is.na(genotype_val), 1, 0),
    AA = if (!is.na(genotype_val) && genotype_val == 2) -1 * AA else AA,
    BB = if (!is.na(genotype_val) && genotype_val == 2) -1 * BB else BB
  ) %>%
  ungroup() %>%
  select(-genotype_val)

peak_summary2 <- filter(peak_summary, Cross_ID == 'GWAS') %>%
  mutate(cross_label = ' GWAS')
peak_summary <- merge(peak_summary2, peak_summary1, 
                      by=intersect(names(peak_summary2), names(peak_summary1)), all=T) %>%
  mutate(chromosome = factor(chromosome, levels = cls$Arm))



##cleaning up missing values manually
peak_summary[peak_summary$chromosome == '7AL' & peak_summary$Cross_ID == "UX1992", c("AA", "BB", "missing_geno")] <- 
  c(peak_summary[peak_summary$chromosome == '7AL' & peak_summary$Cross_ID == "UX1992", c("AA", "BB")]*-1, 0)


# ggplot(peak_summary, aes(x=chromosome, y=cross_label, fill=AA)) +
# # ggplot(peak_summary, aes(x=peak_label, y=cross_label, fill=AA)) +
#   geom_tile(color='white', linewidth=0.2) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
#                        midpoint = 0, na.value = "grey90", name = "AA Effect") +
#   theme_minimal(base_size = 12) +
#   labs(x = "Chromosome", y = "Cross") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text.y = element_text(hjust = 0))


# ggplot(peak_summary, aes(x=chromosome, y=cross_label, fill=AA)) +
#   geom_point(aes(fill = AA, size = LOD), shape=22, colour='gray90') +
#   geom_point(data = subset(peak_summary, missing_geno == 1), aes(shape="Missing parent genotype", color = 'gold', size=2)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
#                        midpoint = 0, na.value = "grey60", name = "AA Effect") +
#   scale_size(range=c(1,8), name='LOD') +
#   scale_shape_manual(
#     name = NULL,
#     values = c("Missing parent genotype" = 8)
#   ) +
#   theme_minimal(base_size = 12) +
#   labs(x = "Chromosome", y = "Cross") +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1),
#         axis.text.y = element_text(hjust = 0))


ggplot(peak_summary, aes(x = chromosome, y = cross_label, fill = AA)) +
  geom_point(aes(fill = AA, size = LOD), shape = 22, colour = 'gray90') +
  
  # Create a factor for the legend and map it to shape
  geom_point(
    data = subset(peak_summary, missing_geno == 1),
    aes(shape = "Missing parent genotype"),
    size = 2
  ) +
  
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, na.value = "grey60", name = "AA Effect"
  ) +
  scale_size(range = c(1, 10), name = 'LOD') +
  scale_shape_manual(
    name = NULL,
    values = c("Missing parent genotype" = 21)
  ) +
  
  theme_minimal(base_size = 12) +
  labs(x = "Chromosome", y = "Cross") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text.y = element_text(hjust = 0),
    legend.box = "vertical",
    legend.position = "right"
  )

ggsave(file="figures/QTL_heatmap_chromArm.png", width=12, height=5)


## generate initial monotonic maps from crosses
# base_map <- read.delim("data/SynOp_RIL906_v1.0_GBS_monotonic.map", header=F)
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
  diff <- 1
  while (diff > 0) {
    r1 <- nrow(lm_mono)
    lm_mono <- lm_mono %>%
      arrange(chrom, cM) %>%
      group_by(chrom) %>%
      mutate(diffp = c(NA, diff(pos)),
             diffc = c(NA, diff(cM))) %>%
      filter(is.na(diffp) | is.na(diffc) | diffp > 0 | diffc > 0) %>%
      select(-c(diffp, diffc)) %>%
      ungroup()
    diff <- r1-nrow(lm_mono)
    print(glue("{r1}, {diff}"))
  }
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
monotonic_consensus <- select(monotonic_consensus, c(chrom, cM, pos))
write.table(monotonic_consensus, "linkage_map/monotonic/consensus_GBS_monotonic.map", quote=F, sep="\t", row.names=F, col.names=F)

for (fam in pedigree$Cross_ID) {
  map <- read.table(glue("linkage_map/monotonic/{fam}_GBS_monotonic.map"), header=F, col.names= c("chrom", "marker", "cM", "pos"))
  for (chr in chroms) {
    if (!(chr %in% map$V1)) {
      add <- monotonic_consensus[monotonic_consensus$chrom == chr,] %>%
        mutate(marker = glue("S{chrom}_{pos}")) %>%
        select(chrom, marker, cM, pos)
      map <- rbind(map, add)
    }
  }
  write.table(map, glue("linkage_map/monotonic/{fam}_GBS_monotonic.map"), quote=F, sep="\t", row.names=F, col.names=F)
  
}
