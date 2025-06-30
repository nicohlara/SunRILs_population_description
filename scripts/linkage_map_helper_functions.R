library("dplyr")
library("qtl")
library("ASMap")
library(stringr)
library(glue)
library(shiny)
library(plotly)
library("ggplot2")
library(pracma)

###notes:
##selection apps should open in browser. When done press 'Done' at bottom of page, then page can be exited and results used
##Phenotype needs to be incorporated BEFORE creating a cross object or it doesn't seem to line up properly

convert_vcf_to_cross <- function(geno, phenotype) {
  pheno <- phenotype %>% filter(genotype %in% row.names(geno))
  pheno <- merge(data.frame(genotype = row.names(geno)), pheno, by='genotype', all=T)
  pheno <- pheno[match(row.names(geno), pheno$genotype),]
  # marker_corr <- apply(geno, 2, function(marker) cor(marker, pheno$Height, use = "pairwise.complete.obs"))
  # flip_markers <- na.omit(names(marker_corr[marker_corr < -0.25]))
  # for (m in flip_markers) {
  #   geno[, m] <- ifelse(is.na(geno[, m]), NA, 2 - geno[, m])  # flips 1 <-> 2
  # }
  chrom <- sapply(strsplit(gsub("S", "", colnames(geno)), "_"), "[", 1)
  geno[geno==0] <- "A"
  geno[geno==1] <- "H"
  geno[geno==2] <- "B"
  geno[is.na(geno)] <- "-"
  geno <- cbind(pheno[1], 1:dim(geno)[1], pheno[2:ncol(pheno)], geno)
  colnames(geno)[2] <- "index"
  geno <- rbind(c(rep("", ncol(pheno)+1),chrom), geno)
  rownames(geno)[1] <- ""
  return(geno)
}  

plot_chrom_curve <- function(cross_map, chrom, zoom=F) {
  plot_frame <- data.frame(cM = cross_map$geno[[chrom]]$map,
                           pos = sapply(strsplit(names(cross_map$geno[[chrom]]$map), "_"), "[", 2))
  par(cex.axis=3, mar=c(5,4,0.5,0.5))
  if (zoom==T) {
    plot(plot_frame$cM, plot_frame$pos, cex=3, pch = 20, xlab="cM", ylab="bp")
  } else {
    plot(plot_frame$cM, plot_frame$pos, ylim=c(1, 1e9), cex=3, pch = 20, xlab="cM", ylab="bp")
  }
}

make_map <- function(map, p_value, min_markers = 10) {
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
  chrnamesDROP<-rownames(subset(summaryMap,summaryMap$n.mar < min_markers))
  markernamesDROP<-markernames(map2,chr=chrnamesDROP)
  map2<-drop.markers(map2,markers=markernamesDROP)
  summaryMap <- summary.map(map2)
  print(summaryMap)
  return(map2)
}

plot_linkage_map <- function(cross_object, family, file_path) {
  summaryMap <- summary.map(cross_object)
  plot_width <- nrow(summaryMap)-1
  plot_depth <- 2
  plot_size <- 4
  text_cex <- 4
  res_val <- 1
  pdf(file=glue("{file_path}{family}_plots.pdf"), width=plot_width*plot_size*res_val, height=plot_depth*plot_size*res_val)
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

extract_group <- function(chrom_name) {
  str_extract(chrom_name, "^[^\\.]+")
}

linkage_group_selector <- function(cross_map, zoom=F) {
  chromosomes <- names(cross_map$geno)
  chrom_groups <- unique(extract_group(chromosomes))
  chrom_by_group <- split(chromosomes, extract_group(chromosomes))
  
  for (chrom in chromosomes) {
    png(filename = file.path("www/plots", paste0(chrom, ".png")), width = 800, height = 600)
    plot_chrom_curve(cross_map, chrom, zoom=zoom)
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


marker_path_selector <- function(df, chr_label = NULL) {
  stopifnot(all(c("cM", "pos") %in% colnames(df)))
  
  env <- new.env()
  env$path_df <- NULL
  
  ui <- fluidPage(
    h4(if (!is.null(chr_label)) paste("Trace a path on chromosome", chr_label) else "Trace a path"),
    plotlyOutput("plot", height = "600px"),
    fluidRow(
      column(4, actionButton("undo", "Undo Last Point")),
      column(4, actionButton("clear", "Clear All")),
      column(4, actionButton("done", "Done"))
    ),
    verbatimTextOutput("status")
  )
  
  server <- function(input, output, session) {
    values <- reactiveValues(path = data.frame(cM = numeric(0), pos = numeric(0)))
    
    observeEvent(event_data("plotly_click", source = "plot_click"), {
      click <- event_data("plotly_click", source = "plot_click")
      if (!is.null(click)) {
        values$path <- rbind(values$path, data.frame(cM = click$x, pos = click$y))
      }
    })
    
    observeEvent(input$undo, {
      if (nrow(values$path) > 0) {
        values$path <- values$path[-nrow(values$path), ]
      }
    })
    
    observeEvent(input$clear, {
      values$path <- data.frame(cM = numeric(0), pos = numeric(0))
    })
    
    output$plot <- renderPlotly({
      plt <- ggplot(df, aes(x = cM, y = pos)) +
        geom_point(alpha = 0.4) +
        geom_path(data = values$path, aes(x = cM, y = pos), color = "blue", linewidth = 1.5) +
        geom_point(data = values$path, aes(x = cM, y = pos), color = "blue", size = 2) +
        theme_minimal()
      
      ggplotly(plt, source = "plot_click") %>%
        layout(dragmode = "click")
    })
    
    output$status <- renderPrint({
      cat("Path points:", nrow(values$path))
    })
    
    observeEvent(input$done, {
      if (nrow(values$path) < 2) {
        env$path_df <- NULL
        stopApp()
      } else {
        path <- values$path
        cM_seq <- seq(min(path$cM), max(path$cM), by = 10)
        interp_pos <- interp1(path$cM, path$pos, cM_seq, method = "linear")
        env$path_df <- data.frame(cM = cM_seq, pos = interp_pos)
        stopApp()
      }
    })
  }
  
  runApp(shinyApp(ui, server), launch.browser = TRUE)
  
  return(env$path_df)
}
