###making final figures for the population manuscript
##Nicolas Lara
##Last edit: 2025-5-19

library(here)
library(tidyverse)
library(gaston)
library(popkin)
library(asreml)
# library(ASRgenomics)
library(yarrr)

setwd(here())

source("C:/Users/nalara/Documents/GitHub/Drone_2023/analysis/palette_functions.R")

###read in complete, imputed, unfiltered vcf file for populations
genotype <- read.vcf("data/SunRILs_prod_filt_imp.vcf.gz", convert.chr=F)
genotype@ped$family <- ifelse(grepl("UX", genotype@ped$id), str_sub(genotype@ped$id, 1, 6), 'Parent')

###read in population information
cross_info <- pedigree <- read.csv("data/cross_info.csv")

###read in BLUPs for traits of interest
blues <- read.delim("data/blues.csv", sep=",") %>%
  mutate(Cross_ID = as.factor(Cross_ID))
##filter
genotype <- select.inds(genotype, family %in% blues$Cross_ID)


###creating palettes and other graphical parameters
light_hue = "white"
dark_hue = "#303D7C"

monochrome_palette <- colorRampPalette(c(light_hue, dark_hue))(100)
index <- c(seq(1,40, by=1), seq(41, 89, by=3), seq(90, 100, by=1))
index <- c(rep(1, 5),
           seq(2, 10, by=1),
           seq(11, 84, by=10),
           seq(90, 100, by=2))
mono_ramp_palette <- monochrome_palette[index]


##ordering families
##set order
parents = c('HILLIARD', "GA06493-13LE6", "NC08-23383")
groups <- data.frame()
for (parent in parents) {
  #groups[[parent]] 
  grp <- cross_info[cross_info$Parent_1 == parent | cross_info$Parent_2 == parent,]
  grp$NAM <- parent
  groups <- rbind(groups, grp)
}
group_order <- c("HILLIARD", "HILLIARD_GA06493-13LE6", "GA06493-13LE6" , "NC08-23383")
groups <- groups %>% group_by(Cross_ID) %>% summarise(NAM = paste(NAM, collapse="_")) %>% 
  mutate(NAM = factor(NAM, levels=group_order)) %>%
  arrange(NAM)

blue_ramp <- generate_ramp(dark_hue, 6)
red_ramp <- generate_ramp('#B71848', 4)
yellow_ramp <- generate_ramp("#9b861b", 4, hue_shift=15)
intermediate <- "#1B731E"
SunRILs_palette <- c(blue_ramp, intermediate, yellow_ramp, red_ramp)
Sun_pal <- list(Hilliard = blue_ramp,
                Hi_GA = intermediate,
                GA = yellow_ramp,
                NC = red_ramp)

plot_palette(SunRILs_palette)
groups$colors <- SunRILs_palette


##table of family contents
genotyped <- genotype@ped %>% group_by(family) %>% count() %>% rename(Cross_ID = family, genotyped = n)
phenotyped <- blues %>% group_by(Cross_ID) %>% count() %>% rename(phenotyped = n)
pop_table <- merge(genotyped, phenotyped, by="Cross_ID" )
pop_table <- merge(cross_info, pop_table, by="Cross_ID")

###Kinship Matrix viewing
genotype_subset <- select.inds(genotype, family != 'Parent')
##create summarized matrix
X <- t(as.matrix(genotype_subset))
kinship <- popkin(X)
k <- aggregate(kinship, list(Family = genotype_subset@ped$family), mean)
rownames(k) <- k[,1]
k <- aggregate(t(k[-1]), list(Family = genotype_subset@ped$family), mean)
rownames(k) <- k[,1]
k <- as.matrix(k[-1])
##set order
k <- k[order(match(rownames(k), groups$Cross_ID), decreasing=T), order(match(colnames(k), groups$Cross_ID), decreasing=T)]
line_plotter <- data.frame(Family = groups$Cross_ID)
join_dataframe <- data.frame(population = c(cross_info$Cross_ID, cross_info$Cross_ID),
                             parent = c(cross_info$Parent_1, cross_info$Parent_2))
pop_index <- data.frame(population = rev(rownames(k)),
                        index = seq(dim(k)[1], 1))
join_dataframe <- merge(join_dataframe, pop_index, by='population' )
parent_order <- c("HILLIARD",
                  "GA00190-7A14",
                  "AGS2000",v 1.3.23
                  "ARGA051160-14LE31",
                  "GA001138-8E36",
                  "GA06493-13LE6",
                  "NC8248-1",
                  "TX12D4896",
                  "SS8641",
                  "MPV57",
                  "LA09264C-P2",
                  "GA05450-EL52",
                  "NC08-23383")
par_index <- data.frame(parent = rev(parent_order),
                        par_index = round(seq(1, max(pop_index$index), by=max(pop_index$index)/13), 0))
join_dataframe <- merge(join_dataframe, par_index, by='parent')
plotter <- select(join_dataframe, c(par_index, index))
plotter <- plotter+2
x_w = 4
x <- c(1, x_w)

res_val=3
png(filename = paste0('figures/kinship_map.png'),  
    width=750*res_val, height=300*res_val, res=72*res_val,
    bg='transparent')
nf <- layout(matrix(c(1,2, 3), 1,3 ,byrow=T),widths=c(x_w+5,ncol(k)+2, x_w-2), heights = c(nrow(k)+2), respect=T)
par(mar=c(1.5,1,.5,0))
plot(x, plotter[1,], ylim=c(0, max(plotter)),
     xlab=NULL, ylab=NULL, axes=F, col='transparent', ann=F, lwd=0)
for (i in rownames(plotter)) {
  lines(x, plotter[i,], lwd=1, col=dark_hue)
}
axis(2, at=par_index$par_index+2, labels=par_index$parent, pos = 1, lty=0, col.axis =dark_hue, las=1, cex.axis=1.5)
par(mar=c(6,0,1,1), cex=.75)
image(k , col=monochrome_palette, axes=F)
axis(side=1, at=seq(0,1, length=15), labels=rownames(k), las=2, cex.axis=1.5, col.axis=dark_hue)
#axis(side=4, at=seq(0,1, length=15), labels=rownames(k), las=2)
image(t(matrix(1:100)), col=monochrome_palette, axes=F)
axis(side=4, at=seq(0, 1, length.out=10), labels=seq(.1,1, length.out=10), col.axis=dark_hue)
mtext('Degree of relatedness', side=4, line=2.5, at=.5,  srt=270, col=dark_hue)
dev.off()

###PCA clustering of populations
M <- as.matrix(genotype)
pcs = prcomp(M)
pc.vars <- data.frame(PC.num = 1:length(pcs$sdev),
                      PC = colnames(pcs$x),
                      var = (pcs$sdev)^2/sum((pcs$sdev)^2))
ggplot(pc.vars,
       aes(x=PC.num,
           y=var)) +
  geom_point() + 
  geom_line() +
  labs(title='Scree plot of PC variances of marker data')
##make classifying group dataframe
line.orig = select(genotype@ped, c(family, id))
data.for.PC.plot = data.frame(pcs$x) %>%
  mutate(id = row.names(pcs$x))
data.for.PC.plot = merge(line.orig, data.for.PC.plot, by='id') %>%
  rename(line = id, group = family) #%>%
data.for.PC.plot <- mutate(data.for.PC.plot, group = factor(group, levels = c(groups$Cross_ID, 'Parent')))
##visualize PCA grouping
HIL <- data.for.PC.plot[data.for.PC.plot$line=="HILLIARD",1:4]
GA <- data.for.PC.plot[data.for.PC.plot$line=="GA06493-13LE6",1:4]
NC08 <- data.for.PC.plot[data.for.PC.plot$line=="NC08-23383",1:4]

ggplot(data.for.PC.plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = group), size=4, alpha = 0.7) +
  #scale_color_viridis_d(option = "H") +
  scale_color_manual(values = c(groups$colors, "black")) +
  theme_minimal() +
  annotate("text", x=HIL[,3] + 2, y=HIL[,4] -3, colour="darkblue", label = "Hilliard") +
  annotate("text", x=GA[,3] + 2, y=GA[,4] - 3, colour="yellow4", label = "GA06493-13LE6") +
  annotate("text", x=NC08[,3] - 3, y=NC08[,4] +6, colour="maroon", label = "NC08-23383") + 
  xlab(paste0(pc.vars[1,'PC'], " ", round(pc.vars[1,'var']*100, 1), "%")) +
  ylab(paste0(pc.vars[2,'PC'], " ", round(pc.vars[2,'var']*100, 1), "%"))
ggsave('figures/pca_plot.png')


##plot phenotypes
parents1 = c('HILLIARD', "GA06493-13LE6", "NC08-23383")
groups <- data.frame()
for (parent in parents1) {
  grp <- cross_info[cross_info$Parent_1 == parent | cross_info$Parent_2 == parent,]
  grp$NAM <- parent
  groups <- rbind(groups, grp)
}
groups <- filter(groups, Cross_ID != 'UX1992')
groups <- rbind(groups, c("UX1992", "GA06493-13LE6", "HILLIARD", "HILLIARD_GA06493-13LE6"))
group_order <- c("HILLIARD", "HILLIARD_GA06493-13LE6", "GA06493-13LE6" , "NC08-23383")
groups <- groups %>% group_by(Cross_ID) %>% summarise(NAM = paste(NAM, collapse="_")) %>% 
  mutate(NAM = factor(NAM, levels=group_order))
cols <- data.frame(NAM = c('GA06493-13LE6', 'HILLIARD', 'NC08-23383', 'HILLIARD_GA06493-13LE6'),
                   cols = c("#B49844", "#313D7D", "#B81848", '#326955'))

groups <- merge(groups, cols, by='NAM') %>%
  arrange(Cross_ID)

res_val=3
png(filename = glue('figures/trait_summaries.png'),  
    width=1700*res_val, height=500*res_val, res=72*res_val,
    bg='transparent')
nf <- layout(matrix(c(1:4), 2,2 ,byrow=T),widths=rep(3.5, 4), heights = rep(1, 4), respect=T)
par(mar=c(2,2,2,2))

for (trait in colnames(blues)[-c(1:2)]) {
  pp <- pirateplot(glue("{trait} ~ Cross_ID"), data=bls, plot=F)
  fam_order <- pp$summary %>% arrange(avg)
  plot_blues <- blues %>%
    filter(Cross_ID != 'Parent') %>%
    mutate(Cross_ID = factor(Cross_ID, levels=fam_order$Cross_ID)) %>%
    arrange(Cross_ID)
  groups <- groups %>% 
    mutate(Cross_ID = factor(Cross_ID, levels=fam_order$Cross_ID)) %>%
    arrange(Cross_ID)
  pirateplot(glue("{trait} ~ Cross_ID"), data=plot_blues,
             bean.f.col = groups$cols,
             main=trait)
}
dev.off()




gc <- merge(groups, cols,by='NAM', all=F) %>%
  filter(Cross_ID %in% pb$Cross_ID) %>%
  arrange(Cross_ID)
pirateplot(Height ~ Cross_ID, data=pb,
           main="Summary of height by cross",
           bar.f.col=gc$cols,
           theme = 2,
           inf.f.o = 0,
           inf.b.o = 0,
           point.o=.2,
           bar.f.o=.5,
           bean.f.o=.4,
           bean.b.o=.2,
           avg.line.o=0,
           point.col='#302D7C')

## to plot as individual NAMs
# UX1992 <- filter(blues, Cross_ID == "UX1992")
# UX1992a <- cbind(UX1992, data.frame(NAM="HILLIARD"))
# UX1992b <- cbind(UX1992, data.frame(NAM="GA06493-13LE6"))
# 
# blues_grp <- filter(blues, Cross_ID != "UX1992") %>%
#   merge(., groups, by="Cross_ID") %>%
#   rbind(., UX1992a, UX1992b)





### calculate heritability
h2 <- function(pheno_file, trait, kinship) {
  pheno_file <- dplyr::select(pheno_file, c(Location, Year, Cross_ID, Entry, !!as.symbol(trait))) %>%
    filter(!is.na(!!as.symbol(trait))) %>%
    droplevels() %>%
    mutate(Env = as.factor(paste0(Location, "_", Year)), 
           Location = as.factor(Location)) %>%
    complete(Env)

  pheno_file <- filter(pheno_file, Entry %in% colnames(kinship))%>%
    mutate(Entry = factor(Entry, levels=colnames(kinship)))

  # ainv <- ainverse(kinship, type="NSD")
  level_num <<- length(unique(pheno_file$Env))
  print(level_num)
  fix.formula <- as.formula(paste0(trait, '~', 1))
  variance <- asreml(fixed = fix.formula,
                     random = ~vm(Entry, kinship) + Location:Year,
                     data=pheno_file,
                     workspace="8gb")
  print(summary(variance)$varcomp)
  herit2 <- asreml::vpredict(variance, herit ~ V2 / (V1 + V2 + V3))
  herit2$trait <- trait
  return(herit2)
}

phenotype <- read.delim("G:/My Drive/Nico_PhD.lnk/data/phenotype/phenotype.csv", sep=",")
locations <- c('Kinston', 'Raleigh')
years <- c('2022', '2023')
effect_variables <- c("Location", "Year", "Cross_ID", "Entry",  "row", "column")
traits <- data.frame(trait = c("WDR", "flowering", "Powdery_mildew", "Height"),
                     type = c("qualitative", "quantitative", "qualitative", "quantitative"))
phenotype <- phenotype %>%
  mutate(Location = as.factor(Location), Cross_ID = as.factor(Cross_ID), Entry=as.factor(Entry), Year = as.factor(Year),
         row = as.factor(row), column = as.factor(column)) %>%
  dplyr::filter(Year %in% years & Cross_ID %in% pedigree$Cross_ID, Location %in% locations) %>%
  dplyr::select(all_of(c(effect_variables, traits$trait))) %>%
  rename(HD = flowering, PM = Powdery_mildew)

g <- select.inds(genotype, id %in% phenotype$Entry)
# kinship <- gaston::GRM(g, autosome.only=F, chunk = 10)
kinship <- G.matrix(M = as.matrix(g), method="VanRaden")$G

h2_df <- data.frame()
traits <- colnames(phenotype)[-c(1:6)]
for (trait in traits) {
  her2 <- h2(filter(phenotype, Entry %in% genotype@ped$id), trait, kinship)
  h2_df <- rbind(h2_df, her2)
  
}
write.csv(h2_df, "outputs/trait_heritability.csv", row.names=F)



## calculate between-trait correlations
corphen <- dplyr::select(phenotype, c(WDR, HD, PM, Height))
cors <- round(cor(corphen, use="pairwise.complete.obs"), 2)
write.csv(cors, "clipboard", row.names=F)
corphen %>%
  summarise_all(list(mean, sd), na.rm=T) %>%
  round(digits=2)


# ###create chromosomal density map
# library(chromoMap)
# chrom_dim <- genotype@snps %>% group_by(chr) %>%
#   summarise(start = 1, end = max(pos)) %>%
#   data.frame()
# chromoMap(list(chrom_dim),list(markers), n_win.factor = 2)#, data_based_color_map = T, data_type = 'numeric')
# 

library(tidyverse)

# Load data
chroms <- paste0( rep(1:7, 3), rep(c("A", "B", "D"), 7))
chrom_lengths <- read.delim("data/chromosome_lengths.tsv")
markers <- genotype@snps %>% select(id, chr, pos) %>% mutate(end = pos+1)
snps <- markers %>%
  select(chr, pos) %>%
  rename(Chr = chr, Position = pos)

# Define bin size
bin_size <- 5e6  # 1 Mb bins

# Create bins per chromosome
bins <- chrom_lengths %>%
  mutate(bin = map2(Chromosome, End, ~ tibble(
    Chr = .x,
    BinStart = seq(1, .y, by = bin_size),
    BinEnd = pmin(seq(1, .y, by = bin_size) + bin_size - 1, .y)
  ))) %>%
  unnest(bin)

# Assign SNPs to bins
snps_binned <- snps %>%
  inner_join(bins, by = "Chr") %>%
  filter(Position >= BinStart & Position <= BinEnd) %>%
  group_by(Chromosome, BinStart, BinEnd) %>%
  summarise(MarkerCount = n(), .groups = 'drop')

# Merge with all bins to include empty ones
heatmap_data <- bins %>%
  left_join(snps_binned, by = c("Chromosome", "BinStart", "BinEnd")) %>%
  mutate(MarkerCount = replace_na(MarkerCount, 0)) %>%
  mutate(Chromosome = factor(Chromosome, levels = rev(unique(Chromosome))))

# Plot
marden <- ggplot(heatmap_data, aes(x = BinStart, y = Chromosome, fill = MarkerCount)) +
  geom_tile(height = 0.9) +
  scale_fill_viridis_c(option = "magma", trans = "log2") +
  labs(x = "Position (bp)", y = "Chromosome", fill = "Marker Count") +
  theme_minimal()
marden
ggsave("figures/marker_density.png", plot = marden, width = 8, height = 3)

###get average distance between markers
average_distances <- genotype@snps %>%
  arrange(chr, pos) %>%
  group_by(chr) %>%
  summarise(
    avg_distance = mean(diff(pos)),
    n_markers = n()
  )
