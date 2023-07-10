suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("glue")
  library("elsasserlib")
  library("apeglm")
  library("wigglescout")
  library("GenomeInfoDb")
  library("ggrepel")
  library("ggpubr")
})

# folders
result_folder = "../results/deseq2_wrapper/"
data_folder = "../data/"
bigwig_folder = "../data/CutNTag_ChIP-Seq/bw/"
cutntag_folder = "../results/cutntag/"

# for subsetting
proms = "../data/CutNTag_ChIP-Seq/bed/ATAC-seq_no_CTCF_no_G4_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed" 
nonproms = "../data/CutNTag_ChIP-Seq/bed/ATAC-seq_no_CTCF_no_G4_WT_without_promoters_missingaszero_skipzero.ordered.regions.bed"

# bigwigs
g4 = c(glue("{bigwig_folder}G4_CnT_E14_NT_R2_RPGC.bigwig"),
       glue("{bigwig_folder}G4_CnT_E14_NT_R1_RPGC.bigwig"))
g4_tmpyp4 = c(glue("{bigwig_folder}G4_CnT_E14_TMPyP4_6h_R2_RPGC.bigwig"),
         glue("{bigwig_folder}G4_CnT_E14_TMPyP4_6h_R1_RPGC.bigwig"))

# DESeq2 wrapper on promoter regions
proms_fc = bw_bed_diff_analysis(
  bwfiles_c1 = g4,
  bwfiles_c2 = g4_tmpyp4,
  bed = proms,
  label_c1 = "g4",
  label_c2 = "g4_tmpyp4", shrink = TRUE
)

proms_fc = as.data.frame(proms_fc)

volc_input = proms_fc %>% mutate(group = case_when(
  log2FoldChange > 5 & padj < 0.0001 ~ "up",
  log2FoldChange < -5 & padj < 0.0001 ~ "down",
  log2FoldChange >= -5 & log2FoldChange <= 5 ~ "unaltered",
  TRUE ~ "non sign."
)) %>%
  mutate(range = paste(seqnames, start, end, sep = "_")) 

labels = volc_input %>% mutate(abs_fc = abs(log2FoldChange)) %>% top_n(., 10, abs_fc) %>% 
  top_n(., -5, padj) %>% pull(range)

volc_input = volc_input %>% mutate(sign_label = ifelse(range %in% labels, range, ""))

labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey", "non sign." = "#9ecae1")
sizes = c("up" = 2, "down" = 2, "unaltered" = 1, "non sign." = 1)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5, "non sign." = 0.5)

# plot
ggplot_volc_proms = volc_input %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             size = group,
             alpha = group,
             color = group)) +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = -log10(0.0001),
             linetype = "dashed") +
  geom_vline(xintercept = c(-5, 5),
             linetype = "dashed") +
  scale_colour_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-30, 30, 10)),  	 
                     limits = c(-30, 30)) +
  scale_y_continuous(breaks = c(seq(0, 40, 10)),  	 
                     limits = c(0, 40)) +
  labs(
    title = "TMPyP4 vs. control G4 signal at promoter regions (DESeq2)",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
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
  ) +
geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc_proms

ggsave(
  glue("{cutntag_folder}wigglescout_TMPyP4_vs_non-tr_DESEq2_volc-proms.pdf"),
  plot = ggplot_volc_proms,
  width = 12,
  height = 8,
  device = "pdf"
)

ggsave(
  glue("{cutntag_folder}wigglescout_TMPyP4_vs_non-tr_DESEq2_volc-proms.png"),
  plot = ggplot_volc_proms,
  width = 12,
  height = 8,
  dpi = 300,
)

# DESeq2 wrapper on non-promoter regions
nonproms_fc = bw_bed_diff_analysis(
  bwfiles_c1 = g4,
  bwfiles_c2 = g4_tmpyp4,
  bed = nonproms,
  label_c1 = "g4",
  label_c2 = "g4_tmpyp4", shrink = TRUE
)

nonproms_fc = as.data.frame(nonproms_fc)

volc_input = nonproms_fc %>% mutate(group = case_when(
  log2FoldChange > 5 & padj < 0.0001 ~ "up",
  log2FoldChange < -5 & padj < 0.0001 ~ "down",
  log2FoldChange >= -5 & log2FoldChange <= 5 ~ "unaltered",
  TRUE ~ "non sign."
)) %>%
  mutate(range = paste(seqnames, start, end, sep = "_")) 

labels = volc_input %>% mutate(abs_fc = abs(log2FoldChange)) %>% top_n(., 10, abs_fc) %>% 
  top_n(., -5, padj) %>% pull(range)

volc_input = volc_input %>% mutate(sign_label = ifelse(range %in% labels, range, ""))

labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey", "non sign." = "#9ecae1")
sizes = c("up" = 2, "down" = 2, "unaltered" = 1, "non sign." = 1)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5, "non sign." = 0.5)

# plot
ggplot_volc_nonproms = volc_input %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             size = group,
             alpha = group,
             color = group)) +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = -log10(0.0001),
             linetype = "dashed") +
  geom_vline(xintercept = c(-5, 5),
             linetype = "dashed") +
  scale_colour_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-30, 30, 10)),  	 
                     limits = c(-30, 30)) +
  scale_y_continuous(breaks = c(seq(0, 40, 10)),  	 
                     limits = c(0, 40)) +
  labs(
    title = "TMPyP4 vs. control G4 signal at non-promoter regions (DESeq2)",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
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
  ) +
  geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc_nonproms

ggsave(
  glue("{cutntag_folder}wigglescout_TMPyP4_vs_non-tr_DESEq2_volc-nonproms.pdf"),
  plot = ggplot_volc_nonproms,
  width = 12,
  height = 8,
  device = "pdf"
)

ggsave(
  glue("{cutntag_folder}wigglescout_TMPyP4_vs_non-tr_DESEq2_volc-nonproms.png"),
  plot = ggplot_volc_nonproms,
  width = 12,
  height = 8,
  dpi = 300,
)

volcanos = ggarrange(ggplot_volc_proms, ggplot_volc_nonproms)
volcanos

ggsave(
  glue("{cutntag_folder}wigglescout_TMPyP4_vs_non-tr_DESEq2_volc.pdf"),
  plot = volcanos,
  width = 17,
  height = 8,
  device = "pdf"
)

ggsave(
  glue("{cutntag_folder}wigglescout_TMPyP4_vs_non-tr_DESEq2_volc.png"),
  plot = volcanos,
  width = 17,
  height = 8,
  dpi = 300,
)

# DESeq2 wrapper: CTCF vs. G4
# bigwigs
g4 = c(glue("{bigwig_folder}G4_0h_AID_5million.norm.bw"),
       glue("{bigwig_folder}G4_6h_AID_5million.norm.bw"))
ctcf = c(glue("{bigwig_folder}CTCF_AID_0h_merge_5million.norm.bw"),
              glue("{bigwig_folder}CTCF_TMPyP4_0h_merge_5million.norm.bw"))


ctcf_fc = bw_bed_diff_analysis(
  bwfiles_c1 = g4,
  bwfiles_c2 = ctcf,
  bed = proms,
  length_factor = 10,
  label_c1 = "g4",
  label_c2 = "g4_ctcf", shrink = TRUE
)

ctcf_fc = as.data.frame(ctcf_fc)

volc_input = ctcf_fc %>% mutate(group = case_when(
  log2FoldChange > 5 & padj < 0.0001 ~ "up",
  log2FoldChange < -5 & padj < 0.0001 ~ "down",
  log2FoldChange >= -5 & log2FoldChange <= 5 ~ "unaltered",
  TRUE ~ "non sign."
)) %>%
  mutate(range = paste(seqnames, start, end, sep = "_")) 

labels = volc_input %>% mutate(abs_fc = abs(log2FoldChange)) %>% top_n(., 10, abs_fc) %>% 
  top_n(., -5, padj) %>% pull(range)

volc_input = volc_input %>% mutate(sign_label = ifelse(range %in% labels, range, ""))

labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey", "non sign." = "#9ecae1")
sizes = c("up" = 2, "down" = 2, "unaltered" = 1, "non sign." = 1)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5, "non sign." = 0.5)

# plot
ggplot_volc_ctcf_proms = volc_input %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             size = group,
             alpha = group,
             color = group)) +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = -log10(0.0001),
             linetype = "dashed") +
  geom_vline(xintercept = c(-5, 5),
             linetype = "dashed") +
  scale_colour_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-30, 30, 10)),  	 
                     limits = c(-30, 30)) +
  scale_y_continuous(breaks = c(seq(0, 40, 10)),  	 
                     limits = c(0, 40)) +
  labs(
    title = "CTCF vs. G4 signal at promoter regions (DESeq2)",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
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
  ) +
  geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc_ctcf_proms

ggsave(
  glue("{cutntag_folder}wigglescout__CTCF_vs_G4_DESEq2_volc-proms.pdf"),
  plot = ggplot_volc_ctcf_proms,
  width = 12,
  height = 8,
  device = "pdf"
)

ggsave(
  glue("{cutntag_folder}wigglescout_CTCF_vs_G4_DESEq2_volc-proms.png"),
  plot = ggplot_volc_ctcf_proms,
  width = 12,
  height = 8,
  dpi = 300,
)

# DESeq2 wrapper on non-promoter regions
nonproms_fc = bw_bed_diff_analysis(
  bwfiles_c1 = g4,
  bwfiles_c2 = ctcf,
  bed = nonproms,
  length_factor = 10,
  label_c1 = "g4",
  label_c2 = "g4_ctcf", shrink = TRUE
)

nonproms_fc = as.data.frame(nonproms_fc)

volc_input = nonproms_fc %>% mutate(group = case_when(
  log2FoldChange > 5 & padj < 0.0001 ~ "up",
  log2FoldChange < -5 & padj < 0.0001 ~ "down",
  log2FoldChange >= -5 & log2FoldChange <= 5 ~ "unaltered",
  TRUE ~ "non sign."
)) %>%
  mutate(range = paste(seqnames, start, end, sep = "_")) 

labels = volc_input %>% mutate(abs_fc = abs(log2FoldChange)) %>% top_n(., 10, abs_fc) %>% 
  top_n(., -5, padj) %>% pull(range)

volc_input = volc_input %>% mutate(sign_label = ifelse(range %in% labels, range, ""))

labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey", "non sign." = "#9ecae1")
sizes = c("up" = 2, "down" = 2, "unaltered" = 1, "non sign." = 1)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5, "non sign." = 0.5)

# plot
ggplot_volc_ctcf_nonproms = volc_input %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             size = group,
             alpha = group,
             color = group)) +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = -log10(0.0001),
             linetype = "dashed") +
  geom_vline(xintercept = c(-5, 5),
             linetype = "dashed") +
  scale_colour_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-30, 30, 10)),  	 
                     limits = c(-30, 30)) +
  scale_y_continuous(breaks = c(seq(0, 40, 10)),  	 
                     limits = c(0, 40)) +
  labs(
    title = "CTCF vs. G4 signal at non-promoter regions (DESeq2)",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
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
  ) +
  geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc_ctcf_nonproms

ggsave(
  glue("{cutntag_folder}wigglescout_CTCF_vs_G4_DESEq2_volc-nonproms.pdf"),
  plot = ggplot_volc_nonproms,
  width = 12,
  height = 8,
  device = "pdf"
)

ggsave(
  glue("{cutntag_folder}wigglescout_CTCF_vs_G4_DESEq2_volc-nonproms.png"),
  plot = ggplot_volc_nonproms,
  width = 12,
  height = 8,
  dpi = 300,
)

volcanos = ggarrange(ggplot_volc_ctcf_proms, ggplot_volc_ctcf_nonproms)
volcanos

ggsave(
  glue("{cutntag_folder}wigglescout_CTCF_vs_G4_DESEq2_volc.pdf"),
  plot = volcanos,
  width = 17,
  height = 8,
  device = "pdf"
)

ggsave(
  glue("{cutntag_folder}wigglescout_CTCF_vs_G4_DESEq2_volc.png"),
  plot = volcanos,
  width = 17,
  height = 8,
  dpi = 300,
)
