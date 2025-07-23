
library(ggplot2)
library(here)
library(dplyr)
library(ggrepel)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(purrr)
setwd(here())

##read in data
chrom_sizes <- read.delim("data/chromosome_lengths.tsv")
mapping <- paste0(rep(1:7, each=3), rep(LETTERS[c(1,2,4)], 7))
asrgwas <- read.delim("outputs/asreml_gwas.tsv") %>%
  dplyr::rename(Position = pos, P.value = p.value) %>%
  mutate(Chromosome = match(chrom, mapping))
asrgwas$model <- "ASRgwas"
blink <- read.delim("outputs/gapit_blink_gwas.tsv")
mlm <- read.delim("outputs/gapit_mlm_gwas.tsv")
farmcpu <- read.delim("outputs/gapit_farmcpu_gwas.tsv")
rrBLUP <- read.delim("outputs/rrBLUP_gwas.tsv") %>%
  dplyr::rename(Position = pos, P.value = p.value) %>%
  mutate(Chromosome = match(chr, mapping))
rrBLUP$model <- "rrBLUP"
# mlm_rerun <- read.delim("outputs/gapit_mlm_gwas_sig_markers_fixed.tsv")



##join together
analysis <- rbind(mlm[,c("Chromosome", "Position", "P.value", "model", "trait")], 
           blink[,c("Chromosome", "Position", "P.value", "model", "trait")],
           farmcpu[,c("Chromosome", "Position", "P.value", "model", "trait")],
           asrgwas[,c("Chromosome", "Position", "P.value", "model", "trait")],
           rrBLUP[,c("Chromosome", "Position", "P.value", "model", "trait")])#,
           # mlm_rerun[,c("Chromosome", "Position", "P.value", "model", "trait")])
analysis$Chromosome <- mapping[mapping=analysis$Chromosome]


##process
traits <- list(flowering = 'HD',
               Powdery_mildew = 'PM')

analysis$trait <- sapply(analysis$trait, function(x) if (x %in% names(traits)) traits[[x]] else x)
analysis <- analysis %>% mutate(trait = factor(trait, levels=c("WDR", "HD", "PM", "Height")))

##subset down to more stringent marker set
bonf_threshold <- (0.05 / 60000) ## total number of SNPs available
analysis <- filter(analysis, P.value <= bonf_threshold)

write.table(analysis, "outputs/combined_GWAS.tsv", quote=F, sep="\t", row.names=F)



# analysis <- filter(analysis, model %in% c("MLM", "MLM_rerun"))
##plot figures
ggplot(analysis, aes(x = Position, y = -log10(P.value), color = model)) +
  geom_point(size=1.5) +
  facet_grid(trait ~ Chromosome, scales="free_y") +
  labs(x = "Position", y = "-log10(P.value)", color = "Model") +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"), axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_x_continuous(breaks=c(0, 4e8, 8e8), limits=c(0, 1e9))
ggsave("figures/gwas_peaks.png", plot = last_plot(), width = 12, height = 4)


## group markers into peaks
assoc <- analysis %>%
  dplyr::rename(chromosome = Chromosome, position = Position, p.value = P.value) %>%
  arrange(trait, chromosome, position) %>%
  mutate(position = position) %>%
  setDT

assign_peak_groups <- function(df, group_distance = 1e6) {
  df <- df %>% arrange(position)
  group_id <- integer(nrow(df))
  current_group <- 1
  group_id[1] <- current_group
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      if (df$position[i] - df$position[i - 1] > group_distance) {
        current_group <- current_group + 1
      }
      group_id[i] <- current_group
    }
  } 
  
  df$peak_group <- group_id
  return(df)
}

# Apply to each (trait, model, chromosome) group
# for (gd  in c(1e7, 3e7, 5e7, 7e7, 9e7)) {
gd <- 3e7
  peaks_clustered <- assoc %>%
    group_by(trait, chromosome) %>%
    group_split() %>%
    map_df(assign_peak_groups, group_distance = gd)
  
  # Summarize peaks
  peak_summary <- data.table(peaks_clustered)[, {
    min_p <- min(p.value)
    peak_pos <- position[which.min(p.value)]
    .(
      # chromosome = data.table::first(chromosome),
      peak_start = min(position),
      peak_end = max(position),
      LOD = -log10(min_p),
      LOD_peak = peak_pos,
      model_count = uniqueN(model),
      nmar = uniqueN(position)
    )
  }, by = .(trait, chromosome, peak_group)]
  # print(gd); print(nrow(peak_summary)); print(table(peak_summary$nmar))
# }

peak_summary <- dplyr::select(peak_summary, -c(peak_group))
write.table(peak_summary, "outputs/peaks.tsv", quote=F, row.names=F, sep="\t")
# peak_summary <- dplyr::filter(data.frame(peak_summary), model_count > 1)
# table(dplyr::filter(data.frame(peak_summary), model_count > 1)$trait)

gwas <- ggplot(data.frame(peak_summary), aes(x = LOD_peak, y = LOD, color = trait)) +
  geom_point(size=1.5) +
  facet_grid(trait ~ chromosome, scales="free_y") +
  labs(x = "Position", y = "-log10(P.value)", color = "Trait") +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"), axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_x_continuous(breaks=c(0, 4e8, 8e8), limits=c(0, 1e9))
ggsave("figures/gwas.png", plot = gwas, width = 12, height = 4)

range_padding <- 10e6

qtl_gr <- GRanges(
  seqnames = paste0("Chr", peak_summary$chromosome),
  ranges = IRanges(start = peak_summary$peak_start-range_padding, end = peak_summary$peak_end + range_padding),
  trait = peak_summary$trait
  )


genes <- import("data/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3")
# genes <- genes[genes$type == "gene"]  # keep only gene entries

# Find overlaps
overlaps <- findOverlaps(qtl_gr, genes)

# Combine QTL and gene info
qtl_gene_matches <- data.table(
  trait = as.character(mcols(qtl_gr)$trait[queryHits(overlaps)]),
  chromosome = as.character(seqnames(qtl_gr)[queryHits(overlaps)]),
  peak_start = start(qtl_gr)[queryHits(overlaps)],
  peak_end = end(qtl_gr)[queryHits(overlaps)],
  gene_id = as.character(mcols(genes)$ID[subjectHits(overlaps)]),
  # gene_chr = as.character(seqnames(genes)[subjectHits(overlaps)]),
  gene_start = start(genes)[subjectHits(overlaps)],
  gene_end = end(genes)[subjectHits(overlaps)]#,
  # gene_name = if ("Name" %in% names(mcols(genes))) as.character(mcols(genes)$Name[subjectHits(overlaps)]) else NA
) %>% unique() %>%
  mutate(chromosome = gsub("Chr", "", chromosome))

output_QTL_summary <- merge(peak_summary, qtl_gene_matches, by=c("trait", "chromosome", "peak_start", "peak_end"), all=T)


##import various sources of gene info
gene_names <- read.delim("data/iwgsc_refseqv2.1_geneID_names.txt") %>%
  dplyr::rename(gene_id = Gene.stable.ID, IWGSC_description = Gene.description) %>%
  select(-"Source..gene.")
TGT_gene_desc <- read.delim("data/TGT_GDTable_20250508230759.csv", sep=",") %>%
  dplyr::rename(gene_id = Gene, TGT_description = Description) %>%
  dplyr::select(-c(Location, Expression))

gene_names <- merge(gene_names, TGT_gene_desc, by=c("gene_id"), all=T)

##merge into named and unnamed gene ranges
named_gene_qtl <- merge(qtl_gene_matches, gene_names, by=c("gene_id"), all=T) %>%
  filter(!(Gene.name == "") & !(trait == "")) %>%
  mutate(peak_start = round(peak_start/1e6, 1), peak_end = round(peak_end/1e6, 1)) %>%
  arrange(trait, chromosome, gene_start)
unnamed_gene_qtl <- merge(qtl_gene_matches, gene_names, by=c("gene_id"), all=T) %>%
  filter((Gene.name == "") & !(trait == ""))
write.table(named_gene_qtl, "outputs/named_qtl_peak_genes.tsv", quote=F, row.names=F, sep="\t")

##generate lollipop plot for publication
qtl_data <- output_QTL_summary %>%
  dplyr::select(trait, chromosome, LOD_peak, LOD) %>%
  dplyr::mutate(chromosome = factor(chromosome, levels = unique(sort(chromosome))),
                trait = factor(trait, levels = c("WDR", "HD", "PM", "Height")), # Optional for grouping color
                LOD_peak = round(LOD_peak/1e6, 1),
                LOD = round(LOD, 1),
                ) %>%
  filter(!(is.na(LOD_peak))) %>% 
  unique()

gene_data <- read.delim("outputs/peaks_annotated.csv", sep="\t") %>% 
  dplyr::filter(!(is.na(gene_start))) %>% 
  dplyr::select(chromosome, gene_start, gene_end, Gene, Status) %>%
  mutate(gene_start = gene_start/1e6, gene_end = gene_end/1e6, text_color = ifelse(Status == "Present", "#000000", "#888888")) %>%
  unique() 
  

##import chromosome sizes
chrom_sizes_mb <- chrom_sizes %>%
  dplyr::mutate(Chromosome = factor(Chromosome, levels = unique(sort(Chromosome))),
         Start = Start / 1e6,
         End = End / 1e6)


# Create the lollipop plot
lolli <- ggplot(qtl_data, aes(x = LOD_peak, y = chromosome)) +
  geom_segment(data = chrom_sizes_mb,
               aes(x = Start, xend = End, y = Chromosome, yend = Chromosome),
               color = "gray90", size = 5) +
  geom_point(data = gene_data, aes(x = gene_start, y = chromosome), color = gene_data$text_color, size = 4, pch = '|') +
  geom_text_repel(data = gene_data,
            aes(x = gene_start, y = chromosome, label = Gene),
            color=gene_data$text_color,
            size = 3, max.overlaps=nrow(qtl_data)) +
  geom_point(aes(color = trait, size = LOD)) +
  scale_size_continuous(range = c(2, 6)) +
  labs(
    x = "Position on Chromosome (Mb)",
    y = "Chromosome",
    size = "LOD Score",
    color = "Trait"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )
lolli
ggsave("figures/gene_lollipop.png", plot = lolli, width = 8, height = 5)


##compare to CIM results
qtl2_cim <- read.delim("outputs/qtl2_cim.tsv") 
qtl2_cim <- filter(qtl2_cim, lod >= -log10(0.25 / 995))
qtl2_cim <- qtl2_cim %>%
  mutate(P.value = 10^-lod, model = "qtl2_cim", Position = as.numeric(gsub("S.._", "", qtl2_cim$Pos))) %>%
  select(chr, Position, P.value, model, trait, cross) %>% 
  dplyr::rename(Chromosome = chr) %>%
  arrange(Chromosome, Position)
##group into peaks
cim_assoc <- qtl2_cim %>%
  arrange(trait, Chromosome, Position) %>% 
  setDT
cim_peaks <- qtl2_cim %>%
  arrange(trait, Chromosome, Position) %>%
  dplyr::rename(position = Position) %>% 
  setDT %>%
  group_by(trait, Chromosome) %>%
  group_split() %>%
  map_df(assign_peak_groups, group_distance = 5e7)
cim_summary <- data.table(cim_peaks)[, {
  min_p <- min(P.value)
  peak_pos <- position[which.min(P.value)]
  .(
    peak_start = min(position),
    peak_end = max(position),
    LOD = -log10(min_p),
    LOD_peak = peak_pos,
    pop_count = uniqueN(cross),
    pops = paste(sort(unique(cross)), collapse=", "),
    nmar = uniqueN(position)
  )
}, by = .(trait, Chromosome, peak_group)]
cim_summary %>% select(-c(peak_group, nmar, pop_count)) %>%
  mutate(peak_start = round(peak_start/1e6,1), peak_end=round(peak_end/1e6, 1), 
         LOD_peak=round(LOD_peak/1e6), LOD=round(LOD/1e6)) %>%
  filter(trait == 'PM') %>%
  data.frame()
