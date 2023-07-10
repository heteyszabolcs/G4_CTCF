# for DiffBind analysis
library("DiffBind")
library("parallel")
library("tidyverse")
library("data.table")

# for volcano
library("ggplot2")
library("ggrepel")
library("ggrastr")
library("glue")
library("ggpubr")

# on remote
# Uppmax working dir
setwd("/proj/snic2020-6-3/SZABOLCS/HiC_G4s/data/CutNTag")

# DiffBind analysis on Uppmax!
g4 <- dba(sampleSheet="diffbind_samplesheet.csv")
g4 <- dba.count(g4)
#readRDS("/proj/snic2020-6-3/SZABOLCS/HiC_G4s/results/cutntag/tmpyp4_cnt-db_count.Rds")
g4 <- dba.normalize(g4)
saveRDS(g4, "/proj/snic2020-6-3/SZABOLCS/HiC_G4s/results/cutntag/DiffBind_TMPyP4_vs_non-tr_norm_output.Rds")
g4 <- dba.contrast(g4, minMembers = 2, contrast = c("Treatment", "TMPyP4", "NT"))
g4 <- dba.analyze(g4, bGreylist = FALSE)
g4 <- dba.report(g4)
result = as.data.frame(g4)
write_tsv(result, "/proj/snic2020-6-3/SZABOLCS/HiC_G4s/results/cutntag/DiffBind_TMPyP4_vs_non-tr_DESEq2_output.tsv")

# on local computer
setwd("./")
result = read.table("../results/cutntag/DiffBind_TMPyP4_vs_non-tr_DESEq2_output.tsv", header = TRUE)
result_folder = "../results/cutntag/"

volc_input = result %>% mutate(group = case_when(
  Fold > 4 & FDR < 0.0001 ~ "up",
  Fold < -4 & FDR < 0.0001 ~ "down",
  Fold >= -4 & Fold <= 4 ~ "unaltered",
  TRUE ~ "non sign."
)) %>%
  mutate(range = paste(seqnames, start, end, sep = "_")) %>% 
  mutate(sign_label = case_when(Fold > 4 & FDR < 0.0001 ~ range,
                                Fold < -4 & FDR < 0.0001 ~ range,
                                Fold >= -4 & Fold <= 4 ~ "",
                                TRUE ~ "non sign."))
labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey", "non sign." = "#9ecae1")
sizes = c("up" = 1, "down" = 1, "unaltered" = 1, "non sign." = 1)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5, "non sign." = 0.5)

# plot
ggplot_volc = volc_input %>%
  ggplot(aes(x = Fold,
             y = -log10(FDR),
             size = group,
             alpha = group,
             color = group)) +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = -log10(0.0001),
             linetype = "dashed") +
  geom_vline(xintercept = c(-4, 4),
             linetype = "dashed") +
  scale_colour_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-15, 15, 5)),  	 
                     limits = c(-15, 15)) +
  labs(
    title = "TMPyP4 vs. control Cut&Tag, DiffBind",
    x = "log2FoldChange",
    y = "-log10 FDR",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, color = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )
  #geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc

ggsave(
  glue("{result_folder}DiffBind_TMPyP4_vs_non-tr_DESEq2_volcano.pdf"),
  plot = ggplot_volc,
  width = 8,
  height = 6,
  device = "pdf"
)

ggsave(
  glue("{result_folder}DiffBind_TMPyP4_vs_non-tr_DESEq2_volcano.png"),
  plot = ggplot_volc,
  width = 8,
  height = 6,
  dpi = 300,
)

# promoters
library("GenomicRanges")
g4_ctcf = fread("../data/CutNTag_ChIP-Seq/bed/promoters_CTCF_G4_common.bed")
g4_ctcf$type = "G4 & CTCF"
g4_ctcf = GRanges(
  seqnames = g4_ctcf$V1,
  ranges = IRanges(
    start = g4_ctcf$V2,
    end = g4_ctcf$V3,
    names = g4_ctcf$type,
  )
)

result$type = "DiffBind"
result_gr = GRanges(
  seqnames = result$seqnames,
  ranges = IRanges(
    start = result$start,
    end = result$end,
    names = result$type,
  )
)

# find overlap with CTCF&G4 promoters
ctcf_g4_ol = findOverlaps(result_gr, g4_ctcf, type = "any", ignore.strand = FALSE)
result_gr = result_gr[queryHits(ctcf_g4_ol)]
result_gr = as_tibble(result_gr)
result_gr = result_gr %>% mutate(range = paste(seqnames, start, end, sep = "_"))
volc_input_ctcf_g4 = volc_input %>% inner_join(., result_gr, by = c("range" = "range"))

# plot
ggplot_volc_ctcf_g4 = volc_input_ctcf_g4 %>%
  ggplot(aes(x = Fold,
             y = -log10(FDR),
             size = group,
             alpha = group,
             color = group)) +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = -log10(0.0001),
             linetype = "dashed") +
  geom_vline(xintercept = c(-4, 4),
             linetype = "dashed") +
  scale_colour_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-15, 15, 5)),  	 
                     limits = c(-15, 15)) +
  labs(
    title = "CTCF & G4 promoters",
    x = "log2FoldChange",
    y = "-log10 FDR",
    fill = " "
  ) +
  guides(alpha = "none", size = "none", color = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )
#geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc_ctcf_g4

g4_only = fread("../data/CutNTag_ChIP-Seq/bed/promoters_G4_no_CTCF.bed")
g4_only$type = "G4 only"
g4_only = GRanges(
  seqnames = g4_only$V1,
  ranges = IRanges(
    start = g4_only$V2,
    end = g4_only$V3,
    names = g4_only$type,
  )
)

result$type = "DiffBind"
result_gr = GRanges(
  seqnames = result$seqnames,
  ranges = IRanges(
    start = result$start,
    end = result$end,
    names = result$type,
  )
)

# find overlap with G4 only promoters
g4_only_ol = findOverlaps(result_gr, g4_only, type = "any", ignore.strand = FALSE)
result_gr = result_gr[queryHits(g4_only_ol)]
result_gr = as_tibble(result_gr)
result_gr = result_gr %>% mutate(range = paste(seqnames, start, end, sep = "_"))
volc_input_g4_only = volc_input %>% inner_join(., result_gr, by = c("range" = "range"))

# plot
ggplot_volc_g4_only = volc_input_g4_only %>%
  ggplot(aes(x = Fold,
             y = -log10(FDR),
             size = group,
             alpha = group,
             color = group)) +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = -log10(0.0001),
             linetype = "dashed") +
  geom_vline(xintercept = c(-4, 4),
             linetype = "dashed") +
  scale_colour_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-15, 15, 5)),  	 
                     limits = c(-15, 15)) +
  labs(
    title = "G4 only promoters",
    x = "log2FoldChange",
    y = "-log10 FDR",
    fill = " "
  ) +
  guides(alpha = "none", size = "none", color = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )
#geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc_g4_only

volcanos = ggarrange(ggplot_volc_g4_only, ggplot_volc_ctcf_g4)
volcanos

ggsave(
  glue("{result_folder}DiffBind_TMPyP4_vs_non-tr_promoter_volcanos.png"),
  plot = volcanos,
  width = 12,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}DiffBind_TMPyP4_vs_non-tr_promoter_volcanos.pdf"),
  plot = volcanos,
  width = 12,
  height = 5,
  device = "pdf"
)


ctcf_only = fread("../data/CutNTag_ChIP-Seq/bed/promoters_CTCF_no_G4.bed")
ctcf_only$type = "CTCF only"
ctcf_only = GRanges(
  seqnames = ctcf_only$V1,
  ranges = IRanges(
    start = ctcf_only$V2,
    end = ctcf_only$V3,
    names = ctcf_only$type,
  )
)

result$type = "DiffBind"
result_gr = GRanges(
  seqnames = result$seqnames,
  ranges = IRanges(
    start = result$start,
    end = result$end,
    names = result$type,
  )
)

# find overlap with CTCF only promoters
ctcf_only_ol = findOverlaps(result_gr, ctcf_only, type = "any", ignore.strand = FALSE)
result_gr = result_gr[queryHits(ctcf_only_ol)]
result_gr = as_tibble(result_gr)
result_gr = result_gr %>% mutate(range = paste(seqnames, start, end, sep = "_"))
volc_input_ctcf_only = volc_input %>% inner_join(., result_gr, by = c("range" = "range"))

# plot
ggplot_volc_ctcf_only = volc_input_ctcf_only %>%
  ggplot(aes(x = Fold,
             y = -log10(FDR),
             size = group,
             alpha = group,
             color = group)) +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = -log10(0.0001),
             linetype = "dashed") +
  geom_vline(xintercept = c(-4, 4),
             linetype = "dashed") +
  scale_colour_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-15, 15, 5)),  	 
                     limits = c(-15, 15)) +
  labs(
    title = "CTCF only promoters",
    x = "log2FoldChange",
    y = "-log10 FDR",
    fill = " "
  ) +
  guides(alpha = "none", size = "none", color = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )
#geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc_ctcf_only
