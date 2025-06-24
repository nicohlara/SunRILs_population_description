
library(here)
library(gaston)
library(popkin)

setwd(here())
pedigree <- read.csv("data/cross_info.csv")


genotype <- read.vcf("data/SunRILs_prod_filt_imp.vcf.gz", convert.chr=F)
geno <- as.matrix(genotype)


###calculate PCA for genotype
pca <- prcomp(geno)
pca$rotation <- -1*pca$rotation
pca$x <- -1*pca$x

##plot PCA results
biplot(pca, scale=0)

##total variance explained by PCA
var <- pcs$sdev^2 / sum(pcs$sdev^2)
head(var)



##plot heatmap of populations
genotype@ped$famid <- ifelse(grepl("UX", genotype@ped$id), str_sub(genotype@ped$id, 1, 6), 'Parent')
genotype_subset <- select.inds(genotype, famid != 'Parent')
parent_lines <- select.inds(genotype, famid == 'Parent')
X <- t(as.matrix(genotype_subset))
parents = c('HILLIARD', "GA06493-13LE6", "NC08-23383")
groups <- list()
for (parent in parents) {
  groups[[parent]] <- pedigree[pedigree$Parent_1 == parent | pedigree$Parent_2 == parent,]
}
groups <- data.frame(Cross_ID = pedigree$Cross_ID, parentage = c('GA06493-13LE6', 'GA06493-13LE6', 'GA06/HIL', 'HILLIARD', 'GA06493-13LE6', 'NC08-23383', 'NC08-23383', 'NC08-23383', 'HILLIARD', 'HILLIARD', 'HILLIARD', 'HILLIARD', 'HILLIARD', 'GA06493-13LE6', 'NC08-23383'))
groups <- merge(data.frame(Cross_ID = genotype@ped$famid), groups, by='Cross_ID')
groups$Entry <- colnames(X)
group_order <- c("HILLIARD", "GA06/HIL", "GA06493-13LE6" , "NC08-23383")
groups <- groups %>% mutate(parentage = factor(parentage, levels=group_order))
groups <- groups[order(groups$parentage), ]

pop_order <- groups$Entry
# applies reordering
X <- X[, order(match(colnames(X), pop_order))]
kinship <- popkin(X, groups$Cross_ID) # calculate kinship from X and optional subpop labels

# ###make tree for phylogeny
# library(ape)
# tree <- rtree( 10 )
# 
# # plot it!
# plot_phylo( tree )
join_dataframe <- data.frame(population = c(pedigree$Cross_ID, pedigree$Cross_ID),
                             parent = c(pedigree$Parent_1, pedigree$Parent_2))
pop_index <- data.frame(entry = rownames(kinship),
                        index = seq(dim(kinship)[1], 1))
pop_fam <- genotype_subset@ped[,c('id', 'famid')] %>%
  mutate(entry = id)
pop_index <- merge(pop_index, pop_fam, by='entry', all=F)
pop_val <- data.frame()
for (i in unique(pop_index$famid)) {
  a <- filter(pop_index, famid == i)
  max <- max(a$index)
  min <- min(a$index)
  pop_val <- rbind(pop_val,
                   data.frame(population = i,
                              pop_loc = round((max+min)/2, 0)))
}
join_dataframe <- merge(join_dataframe, pop_val, by='population' )
parent_order <- c("HILLIARD",
                  "GA00190-7A14",
                  "AGS2000",
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
                        par_index = round(seq(100, max(pop_index$index), by=max(pop_index$index)/13), 0))
join_dataframe <- merge(join_dataframe, par_index, by='parent')
plotter <- select(join_dataframe, c(par_index, pop_loc))
x <- c(1, 1000)





#plot_popkin( inbr_diag(kinship), labs = subpops, labs_las=2, col=c())
full_palette <- colorRampPalette(c("white", "#303D7C"))(100)
index <- c(seq(1,40, by=1), seq(41, 89, by=3), seq(90, 100, by=1))
index <- c(rep(1, 5),
           seq(2, 10, by=1),
           seq(11, 84, by=10),
           seq(90, 100, by=2))
palette <- full_palette[index]

png(width=3500, height=3000, filename = 'figures/kinship_heatmap.png')

par(bg=NA, cex=5, col='#303D7C')

plot_popkin(
  inbr_diag( kinship ),
  labs = cbind(groups$Cross_ID, as.character(groups$parentage) ), # ... labs is now a matrix with levels on columns
  labs_even = c(F, F),      # ... even spacing for first level only
  # labs_line = c(1, 13),             # ... put second level further out
  # labs_cex = c(3, 4),            # ... don't shrink second level
  # labs_sep = c(FALSE, TRUE),       # ... draw lines inside heatmap for second level only
  # labs_las=c(2, 0),                 # push up outer margin ylab "Individuals"
  # mar = 16,                         # increase margins again
  # leg_width = .05, # set legend width relative to plot
  labs_col = '#303D7C',
  col = palette
)
dev.off()

png(width=1500, height=3000, filename = 'figures/kinship_heatmap_pop_labels.png')
par(bg=NA, cex=4, col='#303D7C', las=2, mar=c(1,10,1,1))
plot(x, plotter[1,], ylim=c(0, 2000), 
     xlab=NULL, ylab=NULL, axes=F, col="#303D7C", ann=F, lwd=5)
for (i in rownames(plotter)) {
  lines(x, plotter[i,], lwd=5)
}
axis(2, at=par_index$par_index, labels=par_index$parent, pos = -1, lty=0, col.axis ='#303D7C')
dev.off()