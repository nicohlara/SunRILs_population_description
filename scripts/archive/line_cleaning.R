library(gaston)
library(here)
library(ape)
library(ggtree)
library(plotly)

setwd(here())

vcf_u <- read.vcf("../SunRILs_raw.vcf.gz", convert.chr=F)
vcf_u <- select.snps(vcf_u, chr != "UNKNOWN")
chroms <- unique(vcf_u@snps$chr)
# vcf <- select.snps(vcf, chr %in% sample(chroms)[1:7])
nmarker <- 5000
nind <- 50
totmarkers <- c()
for (k in chroms) {
  markers <- vcf_u@snps[vcf_u@snps$chr == k, 'id']
  totmarkers <- c(totmarkers, markers[ seq(1, length(markers), length.out= nmarker)])
}
vcf1 <- select.snps(vcf_u, id %in% totmarkers)

vcf1@ped$famid <- ifelse(grepl("UX", vcf1@ped$id), sub("\\-.*", "", vcf1@ped$id), "Parent")
totinds <- grep("\\:", vcf1@ped$id, value=T)
for (f in unique(vcf1@ped$famid)) {
  inds <- vcf1@ped[vcf1@ped$famid == f, 'id']
  inds <- grep("\\:", inds, value=T, invert = T)
  totinds <- c(totinds, inds[seq(1, length(inds), length.out=nind)])
}
vcf2 <- select.inds(vcf1, id %in% totinds)

##parent check
grep("NC82", vcf@ped$id, value=T)

NC08_cluster1 <- c(1,2,3,4,6,8,9,10)
NC08_cluster2 <- c(5,7) ##in Hilliard 
GA001138_cluster1 <- c(1,3,2,4,5,7) ##close by is 6
GA001138_cluster2 <- c(8)
TX_cluster1 <- c(2,5,4,7,6,1)
TX_cluster2 <- c(3) #in HILLIARD
NC82_cluster1 <- c(11,12,7,13,8,9) #in UX1989
NC82_cluster2 <- c(10, 3,6,1,11,5,2,4) #in UX1989

##pop fishy
#UX2012-40 clusters with UX2000
#UX2010-142 clusters with UX2012
#UX1989-51:2 clusters with UX2023
#UX1989-51:1 clusters with UX2023 but not close t0 51:2
#UX1997-62:1 clusters in HILLIARD
#UX2031-9977 clusters in HILLIARD
#UX1989-55:2 clusters in HILLIARD
#UX2023-59 clusters with 1989
fishy <- c("UX2012-40", "UX2010-142", "UX1989-51:1", "UX1989-51:2", "UX1997-62:1", 
           "UX2031-9977", "UX1989-55:2", "UX2023-59",
           paste0("NC08-23383", ":", NC08_cluster2),
           paste0("GA001138-8E36", ":",  GA001138_cluster2),
           paste0("TX12D4896", ":",  TX_cluster2))

vcf3 <- select.inds(vcf2, !(id %in% fishy))


fishy2 <- c("UX1993-12","UX2026-73", "UX1992-201:1", "UX1992-201:2", "UX1992-197:1", "UX1992-197:2",
            "NX08-23383:11", "UX2010-112", "UX1989-3:2", "UX2031-999", "UX1994-12",
            "UX2029-21")
vcf4 <- select.inds(vcf1, !(id %in% c(fishy, fishy2, grep("\\:", vcf1@ped$id, value=T))))



# vcf <- vcf4

vcf <- select.inds(vcf2, !(id %in% removal))


chr <- vcf@snps$chr
pos <- vcf@snps$pos
marker_names <- paste0(chr, pos)

snps_num <- as.matrix(vcf)
colnames(snps_num) <- marker_names
N_col <- ncol(snps_num)
N_NAs <- vector(length = N_col)
st_dev <- vector(length = N_col)
for(i in 1:N_col){
  # get the current column
  column_i <- snps_num[, i]
  # get the mean of the current column
  mean_i <- mean(column_i, na.rm = TRUE)
  # get the NAs in the current column
  NAs_i <- which(is.na(column_i))
  # record the number of NAs
  N_NAs[i] <- length(NAs_i)
  # replace the NAs in the current column
  column_i[NAs_i] <- mean_i
  # replace the original column with the
  ## updated columns
  snps_num[, i] <- column_i
  # Record standard deviation after imputation
  st_dev[i] <- sd(column_i)
}

length(st_dev[st_dev == 0])

# Remove invariant columns after imputation
snps_var <- snps_num[, st_dev != 0]


#### Hierarchical clustering ####

# Cluster genotypes by distance
snps_dist<- dist(snps_var, method = "euclidean")
hc <- hclust(snps_dist, method = "ward.D2")
phylo <- as.phylo(hc)

##color clustering
clusters <- cutree(hc, k=15)
famid <- vcf@ped$famid
names(famid) <- vcf@ped$id

tip_metadata <- data.frame(
  label = phylo$tip.label,
  # cluster = as.factor(clusters[phylo$tip.label])
  famid = as.factor(famid[phylo$tip.label])
)

my_colors <- c("orange2", "chartreuse", "green4", "turquoise", "gold", "orangered1", "firebrick1",
               "pink", "purple", "purple3", "royalblue", "paleturquoise1", "navy", "yellow", "maroon")
names(my_colors) <- pedigree$Cross_ID



p <- ggtree(phylo) %<+% tip_metadata +  # %<+% joins metadata
  # geom_tippoint(aes(color = cluster), size = 2) +  # colored tip points
  geom_tippoint(aes(color=famid), size=2) +
  geom_tiplab(aes(label = label), size = 2, hjust = -0.1) +  # labels at end
  theme_tree2() +  # nice axis theme
  scale_color_manual(values=my_colors) +  # use a qualitative color palette
  theme(legend.position = "right")

print(p)
ggsave("figures/full_cladogram.png", scale=1, limitsize=F, width=16, height = 160, dpi=300)

fishy3 <- c("UX2029-8", "UX2029-47", "UX2023-41", "UX2010-68", "UX2010-140", "UX1991-85",
            "UX2031-46", "UX2023-58", "UX1991-66", "UX2026-265", "UX2026-230", "UX2026-206", "UX2026-207",
            "UX1997-9")
fishy4 <- c("UX2029-8", "UX2029-47", "UX2012-18", "UX2031-106", "UX1991-140", "UX1991-85", "UX2031-46",
            "UX2023-58", "UX1991-66", "UX1995-70",  "UX2026-265", "UX1995-31", "UX1997-9",
            "UX1997-62:2", "UX1992-302:1", "UX1992-302:2", "UX1992-20:1", "UX1992-20:2")


vcf5 <- select.inds(vcf1, !(id %in% removal) & (famid == "UX1992" | famid == 'Parent'))
vcf <- vcf5

UX1992_disconnect <- c(319, 294, 8, 336, 273, 303, 9, 216, 203, 56, 8, 16, 56, 211, 316, 319,
                       303, 273, 9, 211, 336, 316)
UX1992_remove <- grep(paste(paste0("UX1992-", c(319, 8, 336, 273, 303, 9, 56, 211, 316 ),":"), collapse="|"), vcf@ped$id, value=T)

removal <- c(fishy, fishy2, fishy3, fishy4, UX1992_remove)

keyfile <- read.delim("data/cladogram_keyfile.tsv")
keyfile <- filter(keyfile, !(FullSampleName %in% removal))
keyfile$FullSampleName <- sub(":.*", "", keyfile$FullSampleName)
keyfile <- dplyr::filter(keyfile, !(Flowcell %in% c("2436449055", "2441621533")))
write.table(keyfile, "data/cleaned_joined_keyfile.tsv", quote = F, row.names=F, sep="\t")
