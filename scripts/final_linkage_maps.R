library("dplyr")
library("qtl")
library("ASMap")
library("gaston")
library(stringr)
library(glue)
library(shiny)
library(plotly)
library(here)

setwd(here())

###notes:
##selection apps should open in browser. When done press 'Done' at bottom of page, then page can be exited and results used
##Phenotype needs to be incorporated BEFORE creating a cross object or it doesn't seem to line up properly


plot_chrom_curve <- function(cross_map, chrom) {
  plot_frame <- data.frame(cM = cross_map$geno[[chrom]]$map,
                           pos = sapply(strsplit(names(cross_map$geno[[chrom]]$map), "_"), "[", 2))
  par(cex.axis=3, mar=c(5,4,0.5,0.5))
  plot(plot_frame$cM, plot_frame$pos, ylim=c(1, 1e9), cex=3, pch = 20, xlab="cM", ylab="bp")
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

convert_vcf_to_cross <- function(gaston_object, phenotype) {
  geno <- as.data.frame(as.matrix(gaston_object), stringAsFactors = F)
  pheno <- phenotype %>% filter(genotype %in% row.names(geno))
  pheno <- merge(data.frame(genotype = row.names(geno)), pheno, by='genotype', all=T)
  pheno <- pheno[match(row.names(geno), pheno$genotype),]
  marker_corr <- apply(geno, 2, function(marker) cor(marker, pheno$Height, use = "pairwise.complete.obs"))
  flip_markers <- na.omit(names(marker_corr[marker_corr < -0.25]))
  for (m in flip_markers) {
    geno[, m] <- ifelse(is.na(geno[, m]), NA, 2 - geno[, m])  # flips 1 <-> 2
  }
  
  geno[geno==0] <- "A"
  geno[geno==1] <- "H"
  geno[geno==2] <- "B"
  geno[is.na(geno)] <- "-"
  geno <- cbind(pheno[2:ncol(pheno)], geno)
  
  geno <- cbind(1:dim(geno)[1], geno)
  colnames(geno)[1] <- "index"
  geno <- rbind(c(rep("", ncol(pheno)),gaston_object@snps$chr), geno)
  rownames(geno)[1] <- ""
  geno <- cbind(rownames(geno), geno)
  colnames(geno)[1] <- 'genotype'
  return(geno)
}  

make_map <- function(map, p_value) {
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

extract_group <- function(chrom_name) {
  str_extract(chrom_name, "^[^\\.]+")
}

linkage_group_selector <- function(cross_map) {
  chromosomes <- names(cross_map$geno)
  chrom_groups <- unique(extract_group(chromosomes))
  chrom_by_group <- split(chromosomes, extract_group(chromosomes))
  
  for (chrom in chromosomes) {
    png(filename = file.path("www/plots", paste0(chrom, ".png")), width = 800, height = 600)
    plot_chrom_curve(cross_map, chrom)
    dev.off()
  }
  
  
  ui <- fluidPage(
    titlePanel("Select Chromosomes"),
    fluidRow(
      uiOutput("chrom_grid"),
      actionButton("done", "Done"),
      verbatimTextOutput("selected")
    )
  )
  
  server <- function(input, output, session) {
    addResourcePath("plots", "www/plots")
    selected <- reactiveVal(character())
    
    output$chrom_grid <- renderUI({
      tagList(
        lapply(names(chrom_by_group), function(group) {
          group_chroms <- chrom_by_group[[group]]
          tags$div(
            tags$h4(group),
            tags$div(
              style = "display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr)); gap: 8px;",
              lapply(group_chroms, function(chrom) {
                uiOutput(outputId = paste0("box_", chrom))
              })
            )
          )
        })
      )
    })
    
    # Dynamically render each chromosome box with colored border
    observe({
      lapply(chromosomes, function(chrom) {
        output[[paste0("box_", chrom)]] <- renderUI({
          is_selected <- chrom %in% selected()
          border_color <- if (is_selected) "blue" else "gray"
          tags$div(
            style = sprintf("border: 3px solid %s; padding: 4px; text-align: center;", border_color),
            actionButton(
              inputId = paste0("btn_", chrom),
              label = tags$img(src = paste0("plots/", chrom, ".png"), width = "100%"),
              style = "padding: 0; border: none; background: none;"
            ),
            tags$p(glue("{chrom}, {length(cross_map$geno[[chrom]]$map)}"), style = "margin: 0; font-size: 12px;"),
            # tags$p(, style = "margin: 0; font-size: 12px;")
          )
        })
      })
    })
    
    # Set up dynamic observers
    observe({
      lapply(chromosomes, function(chrom) {
        observeEvent(input[[paste0("btn_", chrom)]], {
          cur <- selected()
          if (chrom %in% cur) {
            selected(setdiff(cur, chrom))
          } else {
            selected(c(cur, chrom))
          }
        })
      })
    })
    
    output$selected <- renderPrint({
      selected()
    })
    
    observeEvent(input$done, {
      stopApp(selected())
    })
    
  }
  
  result <- runApp(shinyApp(ui, server), launch.browser = TRUE)
  return(result)
}

## refine maps by removing markers
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


### Coerce all vcf files to cross objects
pedigree <- read.delim("data/cross_info.csv", sep=",")
blues <- read.delim("data/blues.csv", sep=",") %>%
  rename(genotype = Entry) %>%
  select(-Cross_ID)

for (fam in pedigree$Cross_ID) {
  print(fam)
  vcf <- read.vcf(glue("linkage_map/biparental_vcf/{fam}_subset.vcf.gz"), convert.chr =FALSE) 
  cross <- convert_vcf_to_cross(vcf, blues)
  write.csv(cross, glue("linkage_map/cross_objects/{fam}_rqtl.csv"), row.names=FALSE)
}


### Actually create the linkage maps
### -------------------------------------------------------------------------###
fam <- "UX1989"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 30000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.5
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
# map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
# map2 <- make_map(map2, p_value = 1e-15)
# keep2 <- linkage_group_selector(map2)
# final_map2 <- subset(map2, keep2)
# final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###


### -------------------------------------------------------------------------###
fam <- "UX1991"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 30000)
clones <- genClones(SunCross1, tol=0.975, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.5
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-14)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- make_map(map2, p_value = 1e-22)
keep2 <- linkage_group_selector(map2)
final_map2 <- subset(map2, keep2)
final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX1992"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 22500)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.65
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-12)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
# ##optional: map2 with remaining chromosome groups
# map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
# map2 <- make_map(map2, p_value = 1e-15)
# keep2 <- linkage_group_selector(map2)
# final_map2 <- subset(map2, keep2)
# final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###


### -------------------------------------------------------------------------###
fam <- "UX1993"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 27500)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.6
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
# map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
# map2 <- make_map(map2, p_value = 1e-15)
# keep2 <- linkage_group_selector(map2)
# final_map2 <- subset(map2, keep2)
# final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX1994"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 20000)
# clones <- genClones(SunCross1, tol=0.95, id="genotype")
# clones$cgd
# SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.55
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
# pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
#                     "bonf", type = "l", cex = 0.25)
# SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
segregation_threshold = 1e-8
SunCross3<-pullCross(SunCross2,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))

SunCross4<-pullCross(SunCross3,type="co.located")

totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- make_map(map2, p_value = 1e-15)
keep2 <- linkage_group_selector(map2)
final_map2 <- subset(map2, keep2)
final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX1995"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 25000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.65
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- make_map(map2, p_value = 1e-15)
keep2 <- linkage_group_selector(map2)
final_map2 <- subset(map2, keep2)
final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX1997"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 30000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.4
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
# map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
# map2 <- make_map(map2, p_value = 1e-15)
# keep2 <- linkage_group_selector(map2)
# final_map2 <- subset(map2, keep2)
# final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX2000"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 25000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.5
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
# map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
# map2 <- make_map(map2, p_value = 1e-15)
# keep2 <- linkage_group_selector(map2)
# final_map2 <- subset(map2, keep2)
# final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX2010"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 25000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.55
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
# pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    # "bonf", type = "l", cex = 0.25)
# SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
segregation_threshold = 1e-8
SunCross4<-pullCross(SunCross3,type="seg.distortion",pars=list(seg.thresh=segregation_threshold))

totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-6)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- make_map(map2, p_value = 1e-15)
keep2 <- linkage_group_selector(map2)
final_map2 <- subset(map2, keep2)
final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)
final_map <- flip.order(final_map, c('5B', '2A'))

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=5, 'hk'))
### -------------------------------------------------------------------------###

##fixed marker missing here
### -------------------------------------------------------------------------###
fam <- "UX2012"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 27500)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.8
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross1$pheno$genotype)*(miss_threshold)]))
totmar(SunCross2)
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- make_map(map2, p_value = 1e-15)
keep2 <- linkage_group_selector(map2)
final_map2 <- subset(map2, keep2)
final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX2013"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 27500)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.75
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross1$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
# map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
# map2 <- make_map(map2, p_value = 1e-15)
# keep2 <- linkage_group_selector(map2)
# final_map2 <- subset(map2, keep2)
# final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX2023"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 30000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.7
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross1$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
# map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
# map2 <- make_map(map2, p_value = 1e-15)
# keep2 <- linkage_group_selector(map2)
# final_map2 <- subset(map2, keep2)
# final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX2026"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 25000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.7
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross1$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- make_map(map2, p_value = 1e-25)
keep2 <- linkage_group_selector(map2)
final_map2 <- subset(map2, keep2)
final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX2029"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 25000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.75
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross1$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
# map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
# map2 <- make_map(map2, p_value = 1e-25)
# keep2 <- linkage_group_selector(map2)
# final_map2 <- subset(map2, keep2)
# final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
fam <- "UX2031"
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 35000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.7
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross1$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
# map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
# map2 <- make_map(map2, p_value = 1e-15)
# keep2 <- linkage_group_selector(map2)
# final_map2 <- subset(map2, keep2)
# final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###








### -------------------------------------------------------------------------###
fam <- ""
SunCross<- read.cross(format="csv",file=glue("linkage_map/cross_objects/{fam}_rqtl.csv"),
                      estimate.map=FALSE, na.strings=c("-","NA"),
                      genotypes=c("A","H","B"), crosstype="riself")
##filter by individuals
pg <- profileGen(SunCross, bychr=F, stat.type=c("miss", "dxo", "xo"), id = 'genotype')
hist(pg$stat$miss)
SunCross1 <- subsetCross(SunCross, ind = pg$stat$miss < 30000)
clones <- genClones(SunCross1, tol=0.95, id="genotype")
clones$cgd
SunCross1 <- fixClones(SunCross1, clones$cgd, consensus = T, id= 'genotype')
##filter by markers
nt.bymar <- ntyped(SunCross1, 'mar')
hist(nt.bymar/length(SunCross1$pheno$genotype))
###Aiming for ~4k markers with >100 markers per chromosome
miss_threshold = 0.5
SunCross2 <- drop.markers(SunCross1, names(nt.bymar[nt.bymar < length(SunCross1$pheno$genotype)*(miss_threshold)]))
SunCross3<-pullCross(SunCross2,type="co.located")
pm <- profileMark(SunCross3, stat.type = c("seg.dist"), crit.val =
                    "bonf", type = "l", cex = 0.25)
SunCross4 <- drop.markers(SunCross3, rownames(pm$marker[pm$marker$crit.val == FALSE,]))
totmar(SunCross); totmar(SunCross1); totmar(SunCross2); totmar(SunCross3); totmar(SunCross4)
##Create linkage maps
map1 <- make_map(SunCross4, p_value = 1e-10)
keep <- linkage_group_selector(map1)
final_map <- subset(map1, keep)
summaryMap <- summary.map(map1)
remove <- rownames(summaryMap)[-grep(paste(sapply(strsplit(keep, "\\."), "[", 1), collapse = "|"), rownames(summaryMap))]
##optional: map2 with remaining chromosome groups
map2 <- subset(SunCross4, chr=unique(head(gsub("\\..*", "", remove), -1)))
map2 <- make_map(map2, p_value = 1e-15)
keep2 <- linkage_group_selector(map2)
final_map2 <- subset(map2, keep2)
final_map <- combineMap(final_map, final_map2, id='genotype')
## polish, clean out rogue markers, rename linkage groups, etc.
names(final_map$geno) <- gsub("\\..*", "", names(final_map$geno))
summary.map(final_map)

markers_to_remove <- c()
for (chr in chrnames(final_map)) {
  message("Selecting markers to remove for chromosome: ", chr)
  selected <- lasso_marker_selector(final_map, chr)
  if (!is.null(selected)) {
    markers_to_remove <- c(markers_to_remove, selected)
  }
}
length(markers_to_remove)

final_map <- drop.markers(final_map, markers_to_remove)
final_map <- quickEst(final_map, map.function = "kosambi")
print(summary.map(final_map)['overall',])
write.table(summary.map(final_map), file = glue("linkage_map/maps/{fam}_stats.csv"), quote=F, row.names=T)
plot_linkage_map(final_map, fam, file_path="linkage_map/maps/")
write.cross(final_map, "csv", filestem=glue("linkage_map/maps/{fam}_linkage_map"))
##test power
scantest <- calc.genoprob(final_map)
plot(scanone(scantest, pheno.col=6, 'hk'))
### -------------------------------------------------------------------------###