library(here)
# library(gaston)
library(dplyr)

setwd(here())
pedigree <- read.csv("data/cross_info.csv")

## OneMap test
library(onemap)
library(glue)
library(grid)

chrom_order <- paste0(rep(1:7, each=3), rep(LETTERS[c(1,2,4)], 7))
cluster_map <- function(group, seq_subset=1) {
  if (length(seq_subset) == 1) {seq_subset<-1:length(group$data.name$CHROM)}
  group_cluster <- data.frame(groups = group$groups,
                              chrom = group$data.name$CHROM[seq_subset])
  a <- table(group_cluster)
  a <- a[rowSums(a) > 0.005*vcf$n.mar,]
  a <- a[apply(a, 1, function(row) sum(row != 0)) <= 3,]
  a <- a[,colSums(a) >5]
  print(a)
  return(a)
}

export_map <- function(group_map, group=T) {
  mars <- group_map$seq.num
  markers <- colnames(group_map$data.name$geno)[mars]
  chrom <- group_map$data.name$CHROM[mars]
  if (group==T) {
    cM <- cumsum(c(0, kosambi(group_map$seq.rf)))
  } else {cM <- NA}
  pos <- group_map$data.name$POS[mars]
  df <- data.frame(marker = markers,
                   chr = chrom,
                   cM_pos = cM,
                   pos = pos)
  return(df)
}

map_curve <- function(exp_map) {
  plot(exp_map$cM_pos, exp_map$pos,
       ylim=c(0, 1e9))
}

matrix_heatmap <- function(matrix_data) {
  color_palette <- colorRampPalette(c("blue", "yellow", "red"))(100)
  zlim <- range(matrix_data)
  
  plot(1, type="n", xlim=c(0, ncol(matrix_data)), ylim=c(0, nrow(matrix_data)), 
       xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  
  # Loop through the matrix and plot a colored rectangle at each position
  for (x in 1:ncol(matrix_data)) {
    for (y in 1:nrow(matrix_data)) {
      # Map matrix value to a color
      value <- matrix_data[y, x]
      color_index <- round((value - zlim[1]) / diff(zlim) * (length(color_palette) - 1)) + 1
      rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col=color_palette[color_index], border=NA)#color_palette[color_index], border=NA)
    }
  }
}

##batch for whole population, round 1
for (fam in cross_info$Cross_ID) {
  dir <- glue("linkage_map/{fam}_map1")
  dir.create(dir)
  # fam = "UX1989"
  vcf <- onemap_read_vcfR(vcf = glue("linkage_map/SunRILs_biparentals/{fam}_subset.vcf.gz"),
                          cross="ri self",
                          parent1 = pedigree[pedigree$Cross_ID==fam, "Parent_1"],
                          parent2 = pedigree[pedigree$Cross_ID==fam, "Parent_2"])
  ##find and filter redundant markers
  bins <- find_bins(vcf, exact=F)
  binned_markers <- create_data_bins(vcf, bins)
  ##test segregation and select non distorted markers
  seg_test <- test_segregation(binned_markers)
  no_dist <- select_segreg(seg_test, distorted = F, numbers = TRUE)
  ##estimate RF
  LOD_thr <- suggest_lod(binned_markers)
  RF <- rf_2pts(binned_markers, max.rf=0.5, LOD=LOD_thr)
  ##make linkage groups
  lg <- make_seq(RF, no_dist)
  LG_group <- group(lg, LOD=LOD_thr)
  
  linkage_groups <-  as.data.frame.matrix(cluster_map(LG_group, no_dist))
  linkage_groups$dom_chrom <- factor(apply(linkage_groups, 1, function(x) {colnames(linkage_groups[which.max(x)])}),
                                     chrom_order)
  linkage_groups <- linkage_groups[order(linkage_groups$dom_chrom), ]
  
  
  
  groups <- row.names(linkage_groups)
  
  res_val=1
  plot_size = 4
  text_cex=4
  plot_depth=3
  plot_width=length(groups)
  pdf(file=glue('{dir}/map_comparison1.pdf'), width=plot_width*plot_size*res_val, height=plot_depth*plot_size*res_val)
  nf <- layout(matrix(c(1:(plot_width*plot_depth)), nrow = plot_depth, ncol = plot_width, byrow=F), 
               heights=matrix(c(rep(plot_size/4,plot_width), rep(plot_size, plot_width*(plot_depth-1))), plot_depth, byrow=T), 
               # widths=rep(c(rep(plot_size, plot_width)), plot_depth),
               widths=matrix(rep(c(rep(plot_size, plot_width)), plot_depth), nrow=plot_depth))
  par(mar=c(0,0,0,0))

  
  
  for (grp_id in groups) {
    LG <- make_seq(LG_group, as.numeric(grp_id))
    
    plot.new()
    text(0.5,0.5, paste0(names(sort(table(LG$data.name$CHROM[LG$seq.num]), decreasing=T)[1]), ", grp ", grp_id), cex=text_cex)
    
    LG_ord <- order_seq(input.seq = LG, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "ug", THRES = LOD_thr,
                        touchdown = T, rm_unlinked=T)
    LG_final <- make_seq(LG_ord, "force")
    grp_map <- export_map(LG_final)
    grp_map$group <- grp_id
    
    map_curve(grp_map)
  
    a <- RF$analysis
    a[upper.tri(a)] <- NA
    a <- a[ LG$seq.num,  LG$seq.num]
    a[upper.tri(a)] <- t(a)[upper.tri(a)]
   
    matrix_heatmap(a)
   
    if (exists("pop_map")) {pop_map <- rbind(pop_map, grp_map)} else {pop_map <- grp_map}
  }
  dev.off()
  
  no_group <- make_seq(RF, LG_group$seq.num[!(LG_group$groups %in% groups)])
  no_group_exp <- export_map(no_group, group=F)
  no_group_exp$group <- NA
  pop_map <- rbind(pop_map, no_group_exp)
  write.table(pop_map, file=glue("{dir}/{fam}_map1"))
  rm(pop_map)
}



# reverse_linkage_group <- function(map_group) {
#   map$cM_pos <- rev(map$cM_pos)
#   return(map)
# }

cluster_map2 <- function(group, seq_subset=1) {
  if (length(seq_subset) == 1) {seq_subset<-1:length(group$data.name$CHROM)}
  group_cluster <- data.frame(groups = group$groups,
                              chrom = group$data.name$CHROM[seq_subset])
  a <- table(group_cluster)
  a <- a[rowSums(a) > 0.002*vcf$n.mar,]
  a <- a[apply(a, 1, function(row) sum(row != 0)) <= 5,]
  a <- a[,colSums(a) >5]
  print(a)
  return(a)
}

## batch for whole population, round 2
for (fam in pedigree$Cross_ID[8:15]) {
  dir <- glue("linkage_map/{fam}_map1")
  map <- read.table(glue("{dir}/{fam}_map1"))
  ##merge groups first
  # for (grp_pair in merging[[fam]]) {map$group <- gsub(grp_pair[2], grp_pair[1], map$group)}
  # for (i in remove_grp[[fam]]) { map <- dplyr::filter(map, !(group %in% i))}
  
  
  
  
  vcf <- onemap_read_vcfR(vcf = glue("linkage_map/SunRILs_biparentals/{fam}_subset.vcf.gz"),
                          cross="ri self",
                          parent1 = pedigree[pedigree$Cross_ID==fam, "Parent_1"],
                          parent2 = pedigree[pedigree$Cross_ID==fam, "Parent_2"])
  
  bins <- find_bins(vcf, exact=F)
  binned_markers <- create_data_bins(vcf, bins)
  ##test segregation and select distorted markers to filter against
  seg_test <- test_segregation(binned_markers)
  no_dist <- select_segreg(seg_test, distorted = T, numbers = TRUE)
  ##estimate RF
  LOD_thr <- suggest_lod(binned_markers) + 2
  RF <- rf_2pts(binned_markers, max.rf=0.5, LOD=LOD_thr)
  # ##make linkage groups
  unmapped_markers <- as.integer(sapply(map[is.na(map$group),"marker"], function(x) grep(x, colnames(RF$analysis))))
  seq <- make_seq(RF, unmapped_markers)
  LG_group <- onemap::group(seq, LOD=LOD_thr)
  
  linkage_groups <-  as.data.frame.matrix(cluster_map2(LG_group, unmapped_markers))
  linkage_groups$dom_chrom <- factor(apply(linkage_groups, 1, function(x) {colnames(linkage_groups[which.max(x)])}),
                                     chrom_order)
  linkage_groups <- linkage_groups[order(linkage_groups$dom_chrom), ]
  groups <- row.names(linkage_groups)
  
  
  res_val=1
  plot_size = 4
  text_cex=4
  plot_depth=3
  plot_width=length(groups)
  pdf(file=glue('{dir}/map_comparison2.pdf'), width=plot_width*plot_size*res_val, height=plot_depth*plot_size*res_val)
  nf <- layout(matrix(c(1:(plot_width*plot_depth)), nrow = plot_depth, ncol = plot_width, byrow=F), 
               heights=matrix(c(rep(plot_size/4,plot_width), rep(plot_size, plot_width*(plot_depth-1))), plot_depth, byrow=T), 
               widths=matrix(rep(c(rep(plot_size, plot_width)), plot_depth), nrow=plot_depth))
  par(mar=c(0,0,0,0))
  for (grp_id in groups) {
    LG <- make_seq(LG_group, as.numeric(grp_id))
    
    plot.new()
    text(0.5,0.5, paste0(names(sort(table(LG$data.name$CHROM[LG$seq.num]), decreasing=T)[1]), ", grp ", grp_id), cex=text_cex)
    
    LG_ord <- order_seq(input.seq = LG, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "ug", THRES = LOD_thr,
                        touchdown = T, rm_unlinked=T)
    LG_final <- make_seq(LG_ord, "force")
    grp_map <- export_map(LG_final)
    grp_map$group <- grp_id
    
    map_curve(grp_map)
    
    a <- RF$analysis
    a[upper.tri(a)] <- NA
    a <- a[ LG$seq.num,  LG$seq.num]
    a[upper.tri(a)] <- t(a)[upper.tri(a)]
    
    matrix_heatmap(a)
    
    if (exists("pop_map")) {pop_map <- rbind(pop_map, grp_map)} else {pop_map <- grp_map}
  }
  dev.off()
  
  # no_group <- make_seq(RF, LG_group$seq.num[!(LG_group$groups %in% groups)])
  # no_group_exp <- export_map(no_group, group=F)
  # no_group_exp$group <- NA
  # pop_map <- rbind(pop_map, no_group_exp)
  
  
  write.table(pop_map, file=glue("{dir}/{fam}_map2"))
  rm(pop_map)
}



## selecting linkage groups from the two batches
round1_remove <- list(UX1989 = list(),
                      UX1991 = list(),
                      UX1992 = list(),
                      UX1993 = list(),
                      UX1994 = list(),
                      UX1995 = list(),
                      UX1997 = list(11),
                      UX2000 = list(),
                      UX2010 = list(),
                      UX2012 = list(),
                      UX2013 = list(),
                      UX2023 = list(),
                      UX2026 = list(),
                      UX2029 = list(4),
                      UX2031 = list())


round2_keep <- list(UX1989 = list(8, 10),
                    UX1991 = list(2, 3, 7, 12, 19),
                    UX1992 = list(5),
                    UX1993 = list(3, 4, 5, 2, 22, 24),
                    UX1994 = list(1, 0),
                    UX1995 = list(3, 0, 1),
                    UX1997 = list(1, 4, 18, 11,3, 19, 22),
                    UX2000 = list(3),
                    UX2010 = list(1, 12, 8, 14),
                    UX2012 = list(2, 0, 5, 7),
                    UX2013 = list(1, 6, 8, 4, 10, 12, 9, 19),
                    UX2023 = list(9, 21, 0, 24),
                    UX2026 = list(7),
                    UX2029 = list(2, 4, 14, 16),
                    UX2031 = list(4, 7, 17))

merge_groups <- list(UX1989 = list(c(1.36, 1.37), c(2.2, 2.6), c(2.18, 1.28), c(1.31, 2.23), c(1.34, 2.25), c(1.38, 2.17)),
                     UX1991 = list(c(1.27, 1.2)),
                     UX1992 = list(c(1.1, 2.3), c(1.15, 1.4), c(1.24, 2.14)),
                     UX1993 = list(c(1.17, 2.6), c(1.3, 2.14), c(1.23, 2.20), c(1.24, 1.5)),
                     UX1994 = list(c(2.3, 1.5), c(1.14, 1.16)),
                     UX1995 = list(c(1.15, 1.12), c(1.3, 2.10, 2.11), c(1.24, 1.25), c(2.15, 1.27)),
                     UX1997 = list(c(2.5, 2.6), c(1.2, 1.15)),
                     UX2000 = list(),
                     UX2010 = list(c(1.12, 2.13), c(2.18, 1.29), c(1.31, 1.11), c(1.17, 1.34), c(1.39, 2.21), c(1.24, 1.25), c(1.39, 1.21)),
                     UX2012 = list(),
                     UX2013 = list(c(2.13, 2.14), c(2.16, 2.17, 2.7), c(2.20, 2.21)),
                     UX2023 = list(c(1.15, 1.16), c(1.4, 1.30)),
                     UX2026 = list(c(2.13, 2.14), c(1.11, 1.25), c(2.16, 1.27)),
                     UX2029 = list(c(1.10, 1.12), c(2.6, 1.13), c(2.11, 2.8)),
                     UX2031 = list(c(2.8, 1.19), c(2.11, 1.24)))


for (fam in pedigree$Cross_ID) {
  dir <- glue("linkage_map/{fam}_map1")
  map1 <- read.table(glue("{dir}/{fam}_map1"))
  map2 <- read.table(glue("{dir}/{fam}_map2"))
  
  ##remove malformed linkage groups from map1
  map1_filt <- dplyr::filter(map1, !(is.na(group)| group %in% round1_remove[[fam]]))
  ##merge partial linkage groups
  map_merge <- map1_filt %>% mutate(group = paste0("1.", group))
  map2_merge <- map2 %>% mutate(group = paste0("2.", group))
  # map_merge <- rbind(map_merge, map2_merge)

  for (grps in merge_groups[[fam]]) {
    m2m <- dplyr::filter(map2_merge, group %in% grps)
    map_merge <- rbind(map_merge, m2m)
    map_merge$group <- gsub(paste0(grps, collapse="|"), paste0(grps[1], "m"), map_merge$group)
  }
  
  
  map2_filt <- dplyr::filter(map2_merge, group %in% paste0("2.", round2_keep[[fam]]))
  map <- rbind(map_merge, map2_filt)
  write.table(map, file=glue("{dir}/{fam}_map3"))
  
}


plot_name=glue('{dir}/map_comparison3.pdf')

plot_chroms_all <- function(map_table, plot_name, 
                            reorder_markers = F, onemap_obj=NULL,
                            res_val=1, plot_size=4, text_cex=4) {
  plot_depth <- ifelse(reorder_markers==T, 3, 2)
  
  linkage_groups <- as.data.frame.matrix(table(map_table[,c("group", "chr")]))
  linkage_groups$dom_chrom <- factor(apply(linkage_groups, 1, function(x) {colnames(linkage_groups[which.max(x)])}),
                                     chrom_order)
  linkage_groups <- linkage_groups[order(linkage_groups$dom_chrom), ]
  groups <- row.names(linkage_groups)
  
  plot_width=length(groups)
  pdf(file=plot_name, width=plot_width*plot_size*res_val, height=plot_depth*plot_size*res_val)
  nf <- layout(matrix(c(1:(plot_width*plot_depth)), nrow = plot_depth, ncol = plot_width, byrow=F), 
               heights=matrix(c(rep(plot_size/4,plot_width), rep(plot_size, plot_width*(plot_depth-1))), plot_depth, byrow=T), 
               widths=matrix(rep(c(rep(plot_size, plot_width)), plot_depth), nrow=plot_depth))
  par(mar=c(2,0,1,0))
  
  
  
  
  for (grp_id in groups) {
    # LG <- make_seq(LG_group, as.numeric(grp_id))
    LG <- dplyr::filter(map_table, group == grp_id)
    plot.new()
    chroms <- sort(table(LG$chr), decreasing=T)
    text(0.5,0.75, paste0(paste0(names(chroms[chroms>sum(chroms)*.1]), collapse=", "), " | grp ", grp_id), cex=text_cex)
    text(0.25,0.25, paste0(chroms[chroms>sum(chroms)*.1], collapse=", "), cex=text_cex-1)
    
    if (reorder_markers == T) {
      LOD_thr <- suggest_lod(onemap_obj)
      LG = make_seq(onemap_obj, as.numeric(grp_id)) ##find markers
      LG_ord <- order_seq(input.seq = LG, n.init = 5,
                          subset.search = "twopt",
                          twopt.alg = "ug", THRES = LOD_thr,
                          touchdown = T, rm_unlinked=T)
      LG_final <- make_seq(LG_ord, "force")
      grp_map <- export_map(LG_final)
      grp_map$group <- grp_id
      map_curve(grp_map)
      RF <- rf_2pts(binned_markers, max.rf=0.5, LOD=LOD_thr)
      a <- RF$analysis
      a[upper.tri(a)] <- NA
      a <- a[ LG$seq.num,  LG$seq.num]
      a[upper.tri(a)] <- t(a)[upper.tri(a)]
      matrix_heatmap(a)
    } else {
      map_curve(LG)
    }

    if (exists("pop_map")) {pop_map <- rbind(pop_map, grp_map)} else {pop_map <- grp_map}
  }
  dev.off()
  return(pop_map)
}




























# calc_LOD_mat <- function(input.seq) {
#   LOD <- matrix(0, length(input.seq$seq.num), length(input.seq$seq.num))
#   k <- matrix(c(rep(input.seq$seq.num[1:(length(input.seq$seq.num))], 
#                     each = length(input.seq$seq.num)), rep(input.seq$seq.num[1:(length(input.seq$seq.num))], 
#                                                            length(input.seq$seq.num))), ncol = 2)
#   k <- k[-which(k[, 1] == k[, 2]), ]
#   k <- t(apply(k, 1, sort))
#   k <- k[-which(duplicated(k)), ]
#   LOD.temp <- input.seq$twopt$analysis[k[, c(1, 2)]]
#   LOD[lower.tri((LOD))] <- LOD.temp
#   LOD[upper.tri(LOD)] <- t(LOD)[upper.tri(LOD)]
#   return(LOD)
# }
# grp <- make_seq(LG_group, grp_id)
# order <- ug(grp, hmm=F)
# rf_graph_table(order, inter=F)
# 
# 
# LG_ord <- order_seq(input.seq = grp, n.init = 5,
#                     subset.search = "twopt",
#                     twopt.alg = "ug", THRES = LOD_thr,
#                     touchdown = T)
# 
# 
# (LG_final <- make_seq(LG_ord, "force"))
# rf_graph_table(LG_final)
# ripple_seq(LG_final, ws = 5 , LOD = LOD_thr)

# 
# LG2 <- make_seq(LG_group, 8)
# LG2_ord <-  order_seq(input.seq = LG2, n.init = 5,
#                       subset.search = "twopt",
#                       twopt.alg = "ug", THRES = LOD_thr,
#                       touchdown = T, rm_unlinked=T)
# (LG2_final <- make_seq(LG2_ord, "force"))
# rf_graph_table(LG2_final)
# # ripple_seq(LG2_final, ws=5)
# 
# 
# 
# 
# LOD <- calc_LOD_mat(make_seq(RF, "all"))
# 
# 
# 
# LOD <- get_mat_rf_in(LG, LOD = T)
# 
# twopt_analysis <- rf_2pts(input.obj = vcf)
# for (chr in unique(vcf$CHROM)) {
#   chr <- "1A"
#   map <- onemap::map(make_seq(twopts_f2, chr))
#   
#   marnames <- colnames(map$data.name$geno)[map$seq.num]
#   df <- data.frame(cM = cumsum(c(0, kosambi(map$seq.rf))),
#                    marker = marnames,
#                    pos = gsub(glue("S{chr}_"), "", marnames))
#   plot(df$cM, df$pos)
# }
# 
# 

# 
# 
# ##try segregation
# f2_test <- test_segregation(vcf)
# plot(f2_test)
# no_dist <- select_segreg(f2_test, distorted = F)
# twopts_f2 <- rf_2pts(input.obj = vcf)
# LOD <- (LOD_sug <- suggest_lod(vcf))
# mark_all_f2 <- make_seq(twopts_f2, "all")
# LGs_f2_simple <- group(mark_all_f2) ##most get put in g1
##clustering didn't respect chrom well
# LGs_upgma <- group_upgma(mark_all_f2, expected.groups = 21, inter = T)
# LGs_f2 <- group(mark_all_f2, LOD = LOD_sug, max.rf = 0.5)
# group_cluster <- data.frame(groups = LGs_f2$groups,
#                             chrom = LGs_f2$data.name$CHROM)
# a <- table(group_cluster)
# a <- a[rowSums(a) > 0.005*vcf$n.mar,]
# a <- a[apply(a, 1, function(row) sum(row != 0)) <= 3,]
# a <- a[,colSums(a) >5]
# 
# ##test linkage groups
# LG1_f2 <- make_seq(LGs_f2, 2)
# LG1_rcd_f2 <- rcd(LG1_f2, hmm = FALSE)
# LG1_rec_f2 <- record(LG1_f2, hmm = FALSE)
# LG1_ug_f2 <- ug(LG1_f2, hmm = FALSE)
# 
# ##plots
# LG1_mds_f2 <- mds_onemap(input.seq = LG1_f2, hmm = F)
# rf_graph_table(LG1_rcd_f2)
# rf_graph_table(LG1_rec_f2)
# rf_graph_table(LG1_ug_f2)
# rf_graph_table(LG1_mds_f2)
# rf_graph_table(LG1_ug_f2, inter = TRUE, html.file = "test.html")
# LG1_f2_ord <- order_seq(input.seq = LG1_ug_f2, n.init = 5,
#                         subset.search = "twopt",
#                         twopt.alg = "rcd", THRES = 3)
# LG1_f2_ord 

# 
# 
# ##by chromosome
# chr1A <- make_seq(twopts_f2, "1A")
# map_1A <- onemap::map(chr1A)
# 
# draw_map(map_1A, names = TRUE, grid = TRUE, cex.mrk = 0.7)
# draw_map2(map_1A)
# 
# map_1A$seq.rf
# map_1A$data.name$POS
# 
# 
# map <- cumsum(c(0, kosambi(map_1A$seq.rf)))
# # map <- cumsum(c(0, map_1A$seq.rf))
# marnames <- colnames(map_1A$data.name$geno)[map_1A$seq.num]
# pos <- gsub("S1A_", "", marnames)
# df <- data.frame(cM = map,
#                  marker = marnames,
#                  pos = pos)
# ##move this to the main pipeline
# tags <- c("-SMTISSUE", "-NWG", "-A+", "-A-")
# genotype@ped$id <- gsub(paste0("\\b(", paste0(tags, collapse="|"), ")\\b"), "", genotype@ped$id)
# ##UX1999 and UX2031 are the same family with switched parents. Renaming now
# genotype@ped$id <- gsub("UX1999", "UX2031", genotype@ped$id)

# ##add family field for split
# genotype@ped$famid <- ifelse(grepl("UX", genotype@ped$id), sub("\\-.*", "", genotype@ped$id), "Parent")
# 
# for (fam in pedigree$Cross_ID) {
#   ##select lines in biparental, plus parents
#   lines <- c(grep(fam, genotype@ped$id, value=T), 
#              pedigree[pedigree$Cross_ID==fam, "Parent_1"],
#              pedigree[pedigree$Cross_ID==fam, "Parent_2"])
#   gs <- select.inds(genotype, id %in% lines)
#   ##refilter for subpop
#   gs <- select.snps(gs, maf > 0.05 &
#                       callrate > 0.8)
#   gs <- select.inds(gs, NAs/ncol(gs) <= 0.5 &
#                       N1/ncol(gs) <= 0.20)
#   ##convert to qtl format
#   geno <- as.data.frame(as.matrix(gs), stringAsFactors = F)
#   geno[geno==0] <- "A"
#   geno[geno==1] <- "H"
#   geno[geno==2] <- "B"
#   geno[is.na(geno)] <- "-"
#   geno <- cbind(1:dim(geno)[1], geno)
#   colnames(geno)[1] <- "index"
#   geno <- rbind(c("",gs@snps$chr), geno)
#   rownames(geno)[1] <- ""
#   geno <- cbind(rownames(geno), geno)
#   colnames(geno)[1] <- 'genotype'
# }