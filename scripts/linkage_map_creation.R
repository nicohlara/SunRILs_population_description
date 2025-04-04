library(here)
library(dplyr)
library(onemap)
library(glue)

setwd(here())
pedigree <- read.csv("data/cross_info.csv")
chrom_order <- paste0(rep(1:7, each=3), rep(LETTERS[c(1,2,4)], 7))



## helper functions
plot_chroms_all <- function(map_table, cluster, plot_name, 
                            reorder_markers = F, onemap_obj=NULL, LOD=NULL,
                            res_val=1, plot_size=4, text_cex=4) {
  plot_depth <- ifelse(reorder_markers==T, 3, 2)

  group_names <- row.names(cluster)
  
  plot_width=length(group_names)
  pdf(file=plot_name, width=plot_width*plot_size*res_val, height=plot_depth*plot_size*res_val)
  nf <- layout(matrix(c(1:(plot_width*plot_depth)), nrow = plot_depth, ncol = plot_width, byrow=F), 
               heights=matrix(c(rep(plot_size/4,plot_width), rep(plot_size, plot_width*(plot_depth-1))), plot_depth, byrow=T), 
               widths=matrix(rep(c(rep(plot_size, plot_width)), plot_depth), nrow=plot_depth))
  par(mar=c(2,0,1,0))
  
  for (grp_id in group_names) {
    print(grp_id)
    # LG <- make_seq(LG_group, as.numeric(grp_id))
    LG <- dplyr::filter(map_table, group == grp_id)
    plot.new()
    chroms <- sort(table(LG$chr), decreasing=T)
    text(0.5,0.75, paste0(paste0(names(chroms[chroms>sum(chroms)*.1]), collapse=", "), " | grp ", grp_id), cex=text_cex)
    text(0.25,0.25, paste0(chroms[chroms>sum(chroms)*.1], collapse=", "), cex=text_cex-1)
    
    if (reorder_markers == T) {
      LG_seq = make_seq(onemap_obj, as.numeric(grp_id)) ##find markers
      LG_ord <- order_seq(input.seq = LG_seq, n.init = 5,
                          subset.search = "twopt",
                          twopt.alg = "ug", THRES = LOD,
                          touchdown = T, rm_unlinked=T)
      LG_final <- make_seq(LG_ord, "force")
      grp_map <- export_map(LG_final)
      grp_map$group <- grp_id
      map_curve(grp_map)
      a <- onemap_obj$twopt$analysis
      # marker_names <- colnames(LG_final$data.name$geno)[LG_final$seq.num]
      a[upper.tri(a)] <- NA
      a <- a[LG$marker, LG$marker]
      a[upper.tri(a)] <- t(a)[upper.tri(a)]
      matrix_heatmap(a)
      if (exists("pop_map")) {pop_map <- rbind(pop_map, grp_map)} else {pop_map <- grp_map}
    } else {
      map_curve(LG)
    }
  }
  dev.off()
  if (exists("pop_map")) {return(pop_map); rm(pop_map)}
}

cluster_map <- function(map, markers, min_marker_percent=0.005, max_chrom=3) {
  marker_limit <- min_marker_percent * markers
  if (marker_limit < 2) {marker_limit=2}
  group_cluster <- data.frame(groups = map$group,
                              chrom = map$chr)
  a <- as.data.frame.matrix(table(group_cluster))
  a <- a[rowSums(a) > marker_limit, , drop=FALSE]
  a <- a[apply(a, 1, function(row) sum(row != 0)) <= max_chrom, , drop=FALSE]
  a <- a[,colSums(a) >5, drop=FALSE]
  a <- a[rowSums(a) > marker_limit, , drop=FALSE]
  a <- a[apply(a, 1, function(row) sum(row != 0)) <= max_chrom, , drop=FALSE]
  a <- a[,colSums(a) >5, drop=FALSE]
  a$dom_chrom <- factor(apply(a, 1, function(x) {colnames(a[which.max(x)])}),
                        chrom_order)
  a <- a[order(a$dom_chrom), ]
  print(a)
  return(a)
}

# export_map <- function(group_map, group=T) {
#   mars <- group_map$seq.num
#   markers <- colnames(group_map$data.name$geno)[mars]
#   chrom <- group_map$data.name$CHROM[mars]
#   if (group==T) {
#     cM <- cumsum(c(0, kosambi(group_map$seq.rf)))
#   } else {cM <- NA}
#   pos <- group_map$data.name$POS[mars]
#   df <- data.frame(marker = markers,
#                    chr = chrom,
#                    cM_pos = cM,
#                    pos = pos)
#   return(df)
# }

map_curve <- function(exp_map) {
  plot(exp_map$cM_pos, exp_map$pos,
       ylim=c(0, 1e9))
}

matrix_heatmap <- function(matrix_data) {
  color_palette <- colorRampPalette(c("#BB0000", "#EEAA66", "#FFFF99"))(100)
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

onemap_loadin <- function(vcf_file, markers=NULL, LOD_adjust=NULL) {
  vcf <- onemap_read_vcfR(vcf = vcf_file,
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
  LOD_thr <- suggest_lod(binned_markers) + LOD_adjust
  RF <- rf_2pts(binned_markers, max.rf=0.5, LOD=LOD_thr)
  ##make linkage groups
  if (length(markers) > 0) {
    lg <- make_seq(RF, as.integer(sapply(markers, function(x) grep(x, colnames(RF$analysis)))))
  } else {lg <- make_seq(RF, no_dist)}
  LG_group <- onemap::group(lg, LOD=LOD_thr)
  selected_markers <- LG_group$marnames
  
  map <- data.frame(marker = selected_markers,
                    chr = sapply(strsplit(selected_markers, "_"), function(x) sub("^S", "", x[1])),
                    pos = sapply(strsplit(selected_markers, "_"), "[", 2),
                    group = LG_group$groups)
  return(list('group'=LG_group, 'map' = map, 'nondistorted_markers'=no_dist, "LOD"=LOD_thr))
}


# make_map <- function(group_obj) {
#   mars <- group_obj$seq.num
#   markers <- colnames(group_obj$data.name$geno)[mars]
#   chrom <- group_obj$data.name$CHROM[mars]
#   pos <- group_obj$data.name$POS[mars]
#   group <- group_obj$groups[mars]
#   df <- data.frame(marker = markers,
#                    chr = chrom,
#                    pos = pos,
#                    group = group)
#   return(df)
# }


### ROUND 1
for (fam in pedigree$Cross_ID) {
  dir <- glue("linkage_map/{fam}_map")
  dir.create(dir)

  LG_group <- onemap_loadin(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), LOD_adjust = 2)
  map <- LG_group$map
  
  cluster <- cluster_map(LG_group$map, LG_group$group$n.mar)

  
  map_exp <- plot_chroms_all(map[map$group %in% row.names(cluster),], 
                             cluster = cluster,
                             plot_name=glue('{dir}/map1.pdf'),
                  reorder_markers=T, onemap_obj=LG_group$group, LOD=LG_group$LOD)

  # for (grp in as.integer(row.names(cluster))) {
    # grp_exp <- export_map(make_seq(LG_group$group, grp), group=F)
    # if (exists("no_group_exp")) {no_group_exp <- rbind(no_group_exp, grp_exp)} else {no_group_exp <- grp_exp}
  # }
  no_group_exp <- dplyr::filter(map, !(marker %in% map_exp$marker))
  no_group_exp$group <- NA
  no_group_exp$cM_pos <- NA
  
  map_exp <- rbind(map_exp, no_group_exp)
  write.table(map_exp, file=glue("{dir}/{fam}_map1"))
  rm(map_exp)
}

## ROUND 1 CLEAN
decisions <- read.delim("linkage_map/map1_decisions.csv", sep=",")

map_clean <- function(decision_row) {
  row_data <- decision_row[-c(1,25)]
  row_data <- row_data[!is.na(row_data)]
  groups <- unlist(strsplit(paste0(row_data, collapse="|"), "|", fixed=TRUE))
  groups <- as.integer(groups[groups != ""])
  fam <- decisions[1,1]
  dir <- glue("linkage_map/{fam}_map")
  map <- read.table(glue("{dir}/{fam}_map1"))
  map$group <- ifelse(map$group %in% groups, map$group, NA)
  return(map)
}

##UX1989
map <- map_clean(decisions[1,])
filter(map, group == 36) %>% arrange(cM_pos) ##appears that data is overlapping
map_curve(filter(map, group == 36 & chr == "6A"))
map$group <- ifelse(map$group == 36 & map$chr == '6A', 37, NA)



 ### ROUND 2--no run for UX1997, UX2013
for (fam in pedigree$Cross_ID[-c(1:11)]) {
  dir <- glue("linkage_map/{fam}_map")
  map1 <- read.table(glue("{dir}/{fam}_map1"))

  ##remove clustered markers
  LG_group <- onemap_loadin(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), map1[is.na(map1$group),"marker"])
  map <- make_map(LG_group$group)
  cluster <- cluster_map(map, LG_group$group$n.mar, min_marker_percent = 0.002, max_chrom = 5)
  
  
  pop_map <- plot_chroms_all(map[map$group %in% row.names(cluster),], plot_name=glue('{dir}/map2.pdf'),
                             reorder_markers=T, onemap_obj=LG_group$group, LOD=LG_group$LOD+2)
  groups <- unique(LG_group$group$groups)
  groups <- groups[!(groups %in% row.names(cluster))]
  for (grp in as.integer(groups)) {
    grp_exp <- export_map(make_seq(LG_group$group, grp), group=F)
    if (exists("no_group_exp")) {no_group_exp <- rbind(no_group_exp, grp_exp)} else {no_group_exp <- grp_exp}
  }
  no_group_exp$group <- NA
  pop_map <- rbind(pop_map, no_group_exp)
  write.table(pop_map, file=glue("{dir}/{fam}_map2"))
  rm(pop_map); rm(no_group_exp)
}


##SELECT LINKAGE GROUPS
## selecting linkage groups from the two batches
map1_keep <- list(UX1989 = list(24, 32, 34, 37),
                      UX1991 = list(7, 10, 15, 18, 24, 25),
                      UX1992 = list(9, 7, 13, 24, 28),
                      UX1993 = list(1, 9),
                      UX1994 = list(),
                      UX1995 = list(),
                      UX1997 = list(),
                      UX2000 = list(),
                      UX2010 = list(),
                      UX2012 = list(),
                      UX2013 = list(),
                      UX2023 = list(),
                      UX2026 = list(),
                      UX2029 = list(),
                      UX2031 = list())
map2_keep <- list(UX1989 = list(6, 10),
                  UX1991 = list(10),
                  UX1992 = list(),
                  UX1993 = list(),
                  UX1994 = list(),
                  UX1995 = list(),
                  UX1997 = list(),
                  UX2000 = list(),
                  UX2010 = list(),
                  UX2012 = list(),
                  UX2013 = list(),
                  UX2023 = list(),
                  UX2026 = list(),
                  UX2029 = list(),
                  UX2031 = list())
group_merge <- list(UX1989 = list(c(1.2, 1.6), c(1.9, 2.5)),
                    UX1991 = list(c(1.28, 1.2)),
                    UX1992 = list(c(1.1,1.6)),
                    UX1993 = list(c(1.7, 2.3), c(2.8, 2.6), c(2.12, 2.10)),
                    UX1994 = list(),
                    UX1995 = list(),
                    UX1997 = list(),
                    UX2000 = list(),
                    UX2010 = list(),
                    UX2012 = list(),
                    UX2013 = list(),
                    UX2023 = list(),
                    UX2026 = list(),
                    UX2029 = list(),
                    UX2031 = list())


for (fam in pedigree$Cross_ID) {
  dir <- glue("linkage_map/{fam}_map")
  map1 <- read.table(glue("{dir}/{fam}_map1"))
  map2 <- read.table(glue("{dir}/{fam}_map2"))
  
  ##select single groups from map1
  map1_select <- dplyr::filter(map1, group %in% map1_keep[[fam]])
  ##select single groups from map2
  map2_select <- dplyr::filter(map2, group %in% map2_keep[[fam]])
  ##merge groups from both
  map_merge <- map1_filt %>% mutate(group = paste0("1.", group))
  map2_merge <- map2 %>% mutate(group = paste0("2.", group))
  for (grps in merge_groups[[fam]]) {
    m2m <- dplyr::filter(map2_merge, group %in% grps)
    map_merge <- rbind(map_merge, m2m)
    map_merge$group <- gsub(paste0(grps, collapse="|"), paste0(grps[1], "m"), map_merge$group)
  }
  ##write output and plot
 

  
  
  map <- rbind(map1_select, map2_select, map_merge)
  write.table(map, file=glue("{dir}/{fam}_map3"))
  
}
