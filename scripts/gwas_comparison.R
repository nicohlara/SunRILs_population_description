
library(ggplot2)
library(here)
library(dplyr)
setwd(here())

##read in data
chrom_sizes <- read.delim("data/chromosome_lengths.tsv")
mapping <- paste0(rep(1:7, each=3), rep(LETTERS[c(1,2,4)], 7))
asrgwas <- read.delim("outputs/asreml_gwas.tsv") %>%
  rename(Position = pos, P.value = p.value) %>%
  mutate(Chromosome = match(chrom, mapping))
asrgwas$model <- "ASRgwas"
blink <- read.delim("outputs/gapit_blink_gwas.tsv")
mlm <- read.delim("outputs/gapit_mlm_gwas.tsv")
farmcpu <- read.delim("outputs/gapit_farmcpu_gwas.tsv")
rrBLUP <- read.delim("outputs/rrBLUP_gwas.tsv") %>%
  rename(Position = pos, P.value = p.value) %>%
  mutate(Chromosome = match(chr, mapping))
rrBLUP$model <- "rrBLUP"
qtl2_cim <- read.delim("outputs/qtl2_cim.tsv") 
qtl2_cim <- qtl2_cim %>%
  mutate(P.value = 10^-lod, model = "qtl2_cim", Position = as.numeric(gsub("S.._", "", qtl2_cim$Pos))) %>%
  select(chr, Position, P.value, model, trait) %>% 
  rename(Chromosome = chr)
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


##create table for presenting
library(dplyr)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

# Prepare your data
assoc <- analysis %>%
  rename(chromosome = Chromosome, position = Position, p.value = P.value) %>%
  arrange(trait, model, chromosome, position) %>%
  setDT

# Group nearby SNPs 
group_distance <- 1e7
assoc[, peak_group := {
  chrom_diff <- chromosome != shift(chromosome, fill = chromosome[1])
  pos_diff <- abs(position - shift(position, fill = position[1])) > group_distance
  is_new_peak <- chrom_diff | pos_diff
  is_new_peak[1] <- TRUE
  cumsum(is_new_peak)
}, by = .(trait, model)]

# Summarize peaks
peak_summary <- assoc[, {
  min_p <- min(p.value)
  peak_pos <- position[which.min(p.value)]
  .(
    chromosome = first(chromosome),
    pos_min = min(position),
    pos_max = max(position),
    LOD = -log10(min_p),
    pos_peak = peak_pos
  )
}, by = .(trait, model, peak_group)]
write.table(peak_summary, "outputs/peaks.tsv", quote=F, row.names=F)

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
  rename(gene_id = Gene.stable.ID) %>%
  select(-"Source..gene.")
named_gene_qtl <- merge(qtl_gene_matches, gene_names, by=c("gene_id")) %>%
  filter(!(Gene.name == "")) %>%
  mutate(peak_chr = gsub('Chr', "", peak_chr)) %>%
  rename(pos_min = peak_start, pos_max = peak_end, chromosome = peak_chr) %>%
  select(-c(Gene.description))

output_QTL_summary <- merge(peak_summary, named_gene_qtl, by=c("trait", "chromosome", "pos_min", "pos_max"), all=T, allow.cartesian=T)
write.table(output_QTL_summary, "outputs/qtl_peaks_gene_names.tsv", quote=F, row.names=F)
