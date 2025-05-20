
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
qtl2_cim <- read.delim("outputs/qtl2_cim.tsv") 
qtl2_cim <- qtl2_cim %>%
  mutate(P.value = 10^-lod, model = "qtl2_cim", Position = as.numeric(gsub("S.._", "", qtl2_cim$Pos))) %>%
  select(chr, Position, P.value, model, trait) %>% 
  dplyr::rename(Chromosome = chr) %>%
  filter(P.value <= 10^-2)
mlm_rerun <- read.delim("outputs/gapit_mlm_gwas_sig_markers_fixed.tsv")



##join together
analysis <- rbind(mlm[,c("Chromosome", "Position", "P.value", "model", "trait")], 
           blink[,c("Chromosome", "Position", "P.value", "model", "trait")],
           farmcpu[,c("Chromosome", "Position", "P.value", "model", "trait")],
           asrgwas[,c("Chromosome", "Position", "P.value", "model", "trait")],
           rrBLUP[,c("Chromosome", "Position", "P.value", "model", "trait")],
           mlm_rerun[,c("Chromosome", "Position", "P.value", "model", "trait")])
analysis$Chromosome <- mapping[mapping=analysis$Chromosome]
analysis <- rbind(analysis, qtl2_cim)


##process
traits <- list(flowering = 'HD',
               Powdery_mildew = 'PM')

analysis$trait <- sapply(analysis$trait, function(x) if (x %in% names(traits)) traits[[x]] else x)
analysis <- analysis %>% mutate(trait = factor(trait, levels=c("WDR", "HD", "PM", "Height")))

##subset down to more stringent marker set
bonf_threshold <- (0.05 / 70000) ## total number of SNPs available
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




# Prepare your data
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
peaks_clustered <- assoc %>%
  group_by(trait, chromosome) %>%
  group_split() %>%
  map_df(assign_peak_groups, group_distance = 1e7)

# Summarize peaks
peak_summary <- data.table(peaks_clustered)[, {
  min_p <- min(p.value)
  peak_pos <- position[which.min(p.value)]
  .(
    chromosome = data.table::first(chromosome),
    pos_min = min(position),
    pos_max = max(position),
    LOD = -log10(min_p),
    pos_peak = peak_pos,
    model_count = uniqueN(model)
  )
}, by = .(trait, chromosome, peak_group)]

write.table(peak_summary, "outputs/peaks.tsv", quote=F, row.names=F, sep="\t")

peak_summary <- dplyr::filter(data.frame(peak_summary), model_count > 1)

ggplot(peak_summary, aes(x = pos_peak, y = LOD, color = trait)) +
  geom_point(size=1.5) +
  facet_grid(trait ~ chromosome, scales="free_y") +
  labs(x = "Position", y = "-log10(P.value)", color = "Model") +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"), axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_x_continuous(breaks=c(0, 4e8, 8e8), limits=c(0, 1e9))


qtl_gr <- GRanges(
  seqnames = paste0("Chr", peak_summary$chromosome),
  ranges = IRanges(start = peak_summary$pos_min, end = peak_summary$pos_max),
  trait = peak_summary$trait
  )


genes <- import("data/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3")
# genes <- genes[genes$type == "gene"]  # keep only gene entries

# Find overlaps
overlaps <- findOverlaps(qtl_gr, genes)

# Combine QTL and gene info
qtl_gene_matches <- data.table(
  trait = as.character(mcols(qtl_gr)$trait[queryHits(overlaps)]),
  peak_chr = as.character(seqnames(qtl_gr)[queryHits(overlaps)]),
  peak_start = start(qtl_gr)[queryHits(overlaps)],
  peak_end = end(qtl_gr)[queryHits(overlaps)],
  gene_id = as.character(mcols(genes)$ID[subjectHits(overlaps)]),
  # gene_chr = as.character(seqnames(genes)[subjectHits(overlaps)]),
  gene_start = start(genes)[subjectHits(overlaps)],
  gene_end = end(genes)[subjectHits(overlaps)]#,
  # gene_name = if ("Name" %in% names(mcols(genes))) as.character(mcols(genes)$Name[subjectHits(overlaps)]) else NA
) %>% unique()

gene_names <- read.delim("data/iwgsc_refseqv2.1_geneID_names.txt") %>%
  dplyr::rename(gene_id = Gene.stable.ID) %>%
  select(-"Source..gene.")
named_gene_qtl <- merge(qtl_gene_matches, gene_names, by=c("gene_id")) %>%
  filter(!(Gene.name == "")) %>%
  mutate(peak_chr = gsub('Chr', "", peak_chr)) %>%
  dplyr::rename(pos_min = peak_start, pos_max = peak_end, chromosome = peak_chr) %>%
  select(-c(Gene.description))

TGT_gene_desc <- read.delim("data/TGT_GDTable_20250508230759.csv", sep=",") %>%
  dplyr::rename(gene_id = Gene, TGT_description = Description)
named_gene_qtl <- merge(named_gene_qtl, TGT_gene_desc, by=c("gene_id"))

output_QTL_summary <- merge(peak_summary, named_gene_qtl, by=c("trait", "chromosome", "pos_min", "pos_max"), all=T, allow.cartesian=T) %>%
  dplyr::select(-c(peak_group, chromosome.1, Expression))
write.table(output_QTL_summary, "outputs/qtl_peaks_gene_names.tsv", quote=F, row.names=F, sep="\t")

##generate lollipop plot for publication
qtl_data <- output_QTL_summary %>%
  dplyr::select(trait, chromosome, pos_peak, LOD, Gene.name, gene_start, gene_end) %>%
  dplyr::rename(peak_position = pos_peak, lod = LOD, gene = Gene.name) %>%
  dplyr::mutate(chromosome = factor(chromosome, levels = unique(sort(chromosome))),
                trait = factor(trait), # Optional for grouping color
                peak_position = round(peak_position/1e6, 1),
                gene_start = round(gene_start/1e6, 1),
                gene_end = round(gene_end/1e6, 1),
                lod = round(lod, 1),
                ) %>%
  unique()

qtl_subset_data <- qtl_data %>% dplyr::filter(!is.na(gene)) %>%
  dplyr::select(chromosome, gene_start, gene_end, gene) %>%
  unique()

##import chromosome sizes
chrom_sizes_mb <- chrom_sizes %>%
  dplyr::mutate(Chromosome = factor(Chromosome, levels = unique(sort(Chromosome))),
         Start = Start / 1e6,
         End = End / 1e6)


# Create the lollipop plot
ggplot(qtl_data, aes(x = peak_position, y = chromosome)) +
  geom_segment(data = chrom_sizes_mb,
               aes(x = Start, xend = End, y = Chromosome, yend = Chromosome),
               color = "gray90", size = 5) +
  geom_point(data = qtl_data, aes(x = gene_start, y = chromosome), color = 'black', size = 4, pch = '|') +
  geom_text_repel(data = qtl_subset_data,
            aes(x = gene_start, y = chromosome, label = gene),
            size = 3, max.overlaps=nrow(qtl_data)) +
  geom_point(aes(color = trait, size = lod)) +
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
