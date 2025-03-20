
library(ggplot2)

mapping <- c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D",
             "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D",
             "7A", "7B", "7D")
asrgwas <- read.delim("outputs/asreml_gwas.tsv") %>%
  rename(Position = pos, P.value = p.value) %>%
  mutate(Chromosome = match(chrom, mapping))
asrgwas$model <- "ASRgwas"
blink <- read.delim("outputs/gapit_blink_gwas.tsv")
mlm <- read.delim("outputs/gapit_mlm_gwas.tsv")
rrBLUP <- read.delim("outputs/rrBLUP_gwas.tsv") %>%
  rename(Position = pos, P.value = p.value) %>%
  mutate(Chromosome = match(chrom, mapping))
rrBLUP$model <- "rrBLUP"

a <- rbind(mlm[,c("Chromosome", "Position", "P.value", "model", "trait")], 
               blink[,c("Chromosome", "Position", "P.value", "model", "trait")], 
                   asrgwas[,c("Chromosome", "Position", "P.value", "model", "trait")],
           rrBLUP[,c("Chromosome", "Position", "P.value", "model", "trait")])
a$Chromosome <- mapping[mapping=a$Chromosome]

# Assuming your dataframe is called df
ggplot(a, aes(x = Position, y = -log10(P.value), color = model)) +
  geom_point(size=3) +
  facet_grid(trait ~ Chromosome, scales="free") +
  labs(x = "Position", y = "-log10(P.value)", color = "Model") +
  theme_bw() # +
  # theme(
  #   strip.text.x = element_text(size = 10),
  #   strip.text.y = element_text(size = 10),
  #   axis.text.x = element_text(angle = 45, hjust = 1),
  #   panel.grid.major = element_line(size = 0.1)
  # )
