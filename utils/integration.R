# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("cowplot")
  library("ggpubr")
  library("ComplexHeatmap")
  library("circlize")
  library("matrixStats")
  library("ggpubr")
})

set.seed(42)

# result folder
result_folder = "../results/hi-c/"

# hi-c data
tad_annot = "../results/hi-c/"
tad_annot = list.files(tad_annot, pattern = 'TAD_annots.tsv', full.names = TRUE)
boundary_annot = "../results/hi-c/"
boundary_annot = list.files(boundary_annot, pattern = '*200kb_annot.tsv', full.names = TRUE)

tad = "../data/Hi-C/bed/"
tad = list.files(tad, pattern = "bed", full.names = TRUE)

# RNA-Seq data
deseq2 = "../results/rna_seq_deseq2/reg_log_norm_counts.tsv"
vst = "../results/rna_seq_deseq2/vst_norm_counts.tsv"
deseq2 = fread(deseq2)
vst = fread(vst)

fc = fread("../results/rna_seq_deseq2/RNA-Seq_treat_vs_contr_DESeq2_fc.tsv")
sign = fc %>% filter(abs(log2FoldChange) > 0.25) %>% pull(gene_name)

# function: join together the two TAD data frames and assigns DESeq2 norm. expression levels to the annotated genes
add_expr_feature = function(tad_annot, tad, expression) {
  bed1 = fread(tad_annot)
  bed2 = fread(tad)
  
  bed1 = bed1 %>% mutate(
    starts = as.character(starts),
    ends = as.character(ends),
    seqnames = as.character(seqnames),
    seqnames = paste0("chr", seqnames)
  )
  bed2 = bed2 %>% mutate(
    start = as.character(start),
    end = as.character(end),
    seqname = as.character(seqnames)
  )
  
  joined = bed1 %>% left_join(., bed2, by = c("starts" = "end", "seqnames" = "seqnames")) %>%
    select(seqnames,
           starts,
           ends,
           gene_symbol,
           TAD_start_log2_ins_score = log2_insulation_score_200000) %>%
    left_join(., bed2, by = c("ends" = "start", "seqnames" = "seqnames")) %>%
    select(
      seqnames,
      starts,
      ends,
      gene_symbol,
      TAD_start_log2_ins_score,
      TAD_end_log2_ins_score = log2_insulation_score_200000
    )
  
  joined_to_expr = joined %>% inner_join(., deseq2, on = "gene_symbol")
  
  return(joined_to_expr)
  
}


# add expression feature to all TAD data tables
# common TADs, non-treated all TADs, non-treated only TADs, TMPyP4 all TADs, TMPyP4 only TADs

## common
common = add_expr_feature(tad_annot = "../results/hi-c/common_200kb_TAD_annots.tsv",
                          tad = "../data/Hi-C/bed/common_fastq_merge_200kb.bed",
                          expression = deseq2)

common_long = pivot_longer(
  common,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)
common_long = common_long %>% drop_na() %>% select(-seqnames,-starts,-ends)

# expression VS. insulation score
common_plot = ggplot(common_long, aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "common TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
common_plot

# highlight genes with higher than abs(0.25) DESeq2 fold change
common_long = common_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "altered", "unaltered"))

common_plot2 = ggplot(common_long, aes(x = NT_1, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = diff.expressed)
  ) +
  scale_fill_manual(values = c("#fc9272", "#f0f0f0"), 
                    labels = c("altered", "unaltered")) +
  geom_point(
    aes(fill = diff.expressed),
    size = 2,
    shape = 21,
    colour = "black",
    data = subset(common_long, diff.expressed == "altered")
  ) +
  labs(
    title = "common TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2InsScore",
    fill = "closest gene"
  ) +
  ylim(-5, 2) +
  scale_x_continuous(breaks = c(seq(0, 20, 2)),       
                     limits = c(0, 20)) +
  guides(color = guide_legend(override.aes = list(fill = 5))) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )
common_plot2

# insulation score vs. log2 fold change

common_fc = common %>% inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  select(gene_symbol, TAD_start_log2_ins_score, TAD_end_log2_ins_score, log2FoldChange, padj)

common_fc = pivot_longer(
  common_fc,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)

common_plot3 = ggplot(common_fc, aes(x = log2FoldChange, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = ins_score_type)
  ) +
  scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("end", "start")) +
  scale_color_manual(values = "black") +
  labs(
    title = "common TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    fill = "TAD position"
  ) +
  scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),       
                     limits = c(-1.5, 1.5)) +
  ylim(-5, 2) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5), reverse = TRUE)) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  stat_cor(method = "pearson", label.x = -1.5, label.y = 2, size = 6)
common_plot3

common_fc = common_fc %>% mutate(diff.expressed = ifelse(padj < 0.05, "1", "0"),
                                 alpha_value = ifelse(padj < 0.05, 1, 0.01))

common_plot4 = ggplot(common_fc %>% arrange(diff.expressed), aes(x = log2FoldChange, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "common TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("padj < 0.05"), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
common_plot4


## non-treated
nt = add_expr_feature(tad_annot = "../results/hi-c/NT_200kb_TAD_annots.tsv",
                      tad = "../data/Hi-C/bed/NT_fastq_merge_bd_200kb.bed",
                      expression = deseq2)

nt_long = pivot_longer(
  nt,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)
nt_long = nt_long %>% drop_na() %>% select(-seqnames,-starts,-ends)

nt_plot = ggplot(nt_long, aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "non-treated TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
nt_plot

nt_long = nt_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "altered", "unaltered"))

nt_plot2 = ggplot(nt_long, aes(x = NT_1, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = diff.expressed)
  ) +
  scale_fill_manual(values = c("#fc9272", "#f0f0f0"), 
                    labels = c("altered", "unaltered")) +
  geom_point(
    aes(fill = diff.expressed),
    size = 2,
    shape = 21,
    colour = "black",
    data = subset(common_long, diff.expressed == "altered")
  ) +
  labs(
    title = "non-treated TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2InsScore",
    fill = "closest gene"
  ) +
  ylim(-5, 2) +
  scale_x_continuous(breaks = c(seq(0, 20, 2)),       
                     limits = c(0, 20)) +
  guides(color = guide_legend(override.aes = list(fill = 5))) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )
nt_plot2

# insulation score vs. log2 fold change

nt_fc = nt %>% inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  select(gene_symbol, TAD_start_log2_ins_score, TAD_end_log2_ins_score, log2FoldChange, padj)

nt_fc = pivot_longer(
  nt_fc,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)

nt_plot3 = ggplot(nt_fc, aes(x = log2FoldChange, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = ins_score_type)
  ) +
  scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("end", "start")) +
  labs(
    title = "non-treated TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    fill = "TAD position"
  ) +
  scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),       
                     limits = c(-1.5, 1.5)) +
  ylim(-5, 2) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5), reverse = TRUE)) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  stat_cor(method = "pearson", label.x = -1.5, label.y = 2, size = 6)
nt_plot3

nt_fc = nt_fc %>% mutate(diff.expressed = ifelse(padj < 0.05, "1", "0"),
                                 alpha_value = ifelse(padj < 0.05, 1, 0.01))

nt_plot4 = ggplot(nt_fc %>% arrange(diff.expressed), aes(x = log2FoldChange, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "non-treated TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("padj < 0.05"), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
nt_plot4

## non-treated only
nt_only = add_expr_feature(tad_annot = "../results/hi-c/NT_only_200kb_TAD_annots.tsv",
                      tad = "../data/Hi-C/bed/NT_only_fastq_merge_200kb.bed",
                      expression = deseq2)

nt_only_long = pivot_longer(
  nt_only,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)
nt_only_long = nt_only_long %>% drop_na() %>% select(-seqnames,-starts,-ends)

nt_only_plot = ggplot(nt_only_long, aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "non-treated only TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
nt_only_plot

nt_only_long = nt_only_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "altered", "unaltered"))

nt_only_plot2 = ggplot(nt_only_long, aes(x = NT_1, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = diff.expressed)
  ) +
  scale_fill_manual(values = c("#fc9272", "#f0f0f0"), 
                    labels = c("altered", "unaltered")) +
  geom_point(
    aes(fill = diff.expressed),
    size = 2,
    shape = 21,
    colour = "black",
    data = subset(common_long, diff.expressed == "altered")
  ) +
  labs(
    title = "non-treated only TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2InsScore",
    fill = "closest gene"
  ) +
  ylim(-5, 2) +
  scale_x_continuous(breaks = c(seq(0, 20, 2)),       
                     limits = c(0, 20)) +
  guides(color = guide_legend(override.aes = list(fill = 5))) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )
nt_only_plot2

# insulation score vs. log2 fold change
nt_only_fc = nt_only %>% inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  select(gene_symbol, TAD_start_log2_ins_score, TAD_end_log2_ins_score, log2FoldChange, padj)

nt_only_fc = pivot_longer(
  nt_only_fc,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)

nt_only_plot3 = ggplot(nt_only_fc, aes(x = log2FoldChange, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = ins_score_type)
  ) +
  scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("end", "start")) +
  labs(
    title = "non-treated only TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    fill = "TAD position"
  ) +
  scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),       
                     limits = c(-1.5, 1.5)) +
  ylim(-5, 2) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5), reverse = TRUE)) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  stat_cor(method = "pearson", label.x = -1.5, label.y = 2, size = 6)
nt_only_plot3

nt_only_fc = nt_only_fc %>% mutate(diff.expressed = ifelse(padj < 0.05, "1", "0"),
                                 alpha_value = ifelse(padj < 0.05, 1, 0.01))

nt_only_plot4 = ggplot(nt_only_fc %>% arrange(diff.expressed), aes(x = log2FoldChange, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "non-treated TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("padj < 0.05"), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
nt_only_plot4

## TMPyP4
TMPyP4 = add_expr_feature(tad_annot = "../results/hi-c/TMPyP4_200kb_TAD_annots.tsv",
                           tad = "../data/Hi-C/bed/TMPyP4_fastq_merge_bd_200kb.bed",
                           expression = deseq2)

TMPyP4_long = pivot_longer(
  TMPyP4,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)
TMPyP4_long = TMPyP4_long %>% drop_na() %>% select(-seqnames,-starts,-ends)

TMPyP4_plot = ggplot(TMPyP4_long, aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "TMPyP4-treated TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
TMPyP4_plot

TMPyP4_long = TMPyP4_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "altered", "unaltered"))

TMPyP4_plot2 = ggplot(TMPyP4_long, aes(x = NT_1, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = diff.expressed)
  ) +
  scale_fill_manual(values = c("#fc9272", "#f0f0f0"), 
                    labels = c("altered", "unaltered")) +
  geom_point(
    aes(fill = diff.expressed),
    size = 2,
    shape = 21,
    colour = "black",
    data = subset(common_long, diff.expressed == "altered")
  ) +
  labs(
    title = "TMPyP4 TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2InsScore",
    fill = "closest gene"
  ) +
  ylim(-5, 2) +
  scale_x_continuous(breaks = c(seq(0, 20, 2)),       
                     limits = c(0, 20)) +
  guides(color = guide_legend(override.aes = list(fill = 5))) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )
TMPyP4_plot2

# insulation score vs. log2 fold change
TMPyP4_fc = TMPyP4 %>% inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  select(gene_symbol, TAD_start_log2_ins_score, TAD_end_log2_ins_score, log2FoldChange, padj)

TMPyP4_fc = pivot_longer(
  TMPyP4_fc,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)

TMPyP4_plot3 = ggplot(TMPyP4_fc, aes(x = log2FoldChange, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = ins_score_type)
  ) +
  scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("end", "start")) +
  labs(
    title = "TMPyP4 TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    fill = "TAD position"
  ) +
  scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),       
                     limits = c(-1.5, 1.5)) +
  ylim(-5, 2) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5), reverse = TRUE)) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  stat_cor(method = "pearson", label.x = -1.5, label.y = 2, size = 6)
TMPyP4_plot3

TMPyP4_fc = TMPyP4_fc %>% mutate(diff.expressed = ifelse(padj < 0.05, "1", "0"),
                                   alpha_value = ifelse(padj < 0.05, 1, 0.01))

TMPyP4_plot4 = ggplot(TMPyP4_fc %>% arrange(diff.expressed), aes(x = log2FoldChange, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "TMPyP4-treated TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("padj < 0.05"), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
TMPyP4_plot4

## TMPyP4 only
TMPyP4_only = add_expr_feature(tad_annot = "../results/hi-c/TMPyP4_only_200kb_TAD_annots.tsv",
                          tad = "../data/Hi-C/bed/TMPyP4_only_fastq_merge_200kb.bed",
                          expression = deseq2)

TMPyP4_only_long = pivot_longer(
  TMPyP4_only,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)
TMPyP4_only_long = TMPyP4_only_long %>% drop_na() %>% select(-seqnames,-starts,-ends)

TMPyP4_only_plot = ggplot(TMPyP4_only_long, aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "TMPyP4-treated only TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
TMPyP4_only_plot

TMPyP4_only_long = TMPyP4_only_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "altered", "unaltered")) 

TMPyP4_only_plot2 = ggplot(TMPyP4_only_long, aes(x = NT_1, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = diff.expressed)
  ) +
  scale_fill_manual(values = c("#fc9272", "#f0f0f0"), 
                    labels = c("altered", "unaltered")) +
  geom_point(
    aes(fill = diff.expressed),
    size = 2,
    shape = 21,
    colour = "black",
    data = subset(common_long, diff.expressed == "altered")
  ) +
  labs(
    title = "TMPyP4 only TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2InsScore",
    fill = "closest gene"
  ) +
  ylim(-5, 2) +
  scale_x_continuous(breaks = c(seq(0, 20, 2)),       
                     limits = c(0, 20)) +
  guides(color = guide_legend(override.aes = list(fill = 5))) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )
TMPyP4_only_plot2

# insulation score vs. log2 fold change
TMPyP4_only_fc = TMPyP4_only %>% inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  select(gene_symbol, TAD_start_log2_ins_score, TAD_end_log2_ins_score, log2FoldChange, padj)

TMPyP4_only_fc = pivot_longer(
  TMPyP4_only_fc,
  cols = c("TAD_start_log2_ins_score", "TAD_end_log2_ins_score"),
  names_to = "ins_score_type",
  values_to = "log2_ins_score"
)

TMPyP4_only_plot3 = ggplot(TMPyP4_only_fc, aes(x = log2FoldChange, y = log2_ins_score)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = ins_score_type)
  ) +
  scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("end", "start")) +
  labs(
    title = "TMPyP4 only TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    fill = "TAD position"
  ) +
  scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),       
                     limits = c(-1.5, 1.5)) +
  ylim(-5, 2) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5), reverse = TRUE)) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  stat_cor(method = "pearson", label.x = -1.5, label.y = 2, size = 6)
TMPyP4_only_plot3

TMPyP4_only_fc = TMPyP4_only_fc %>% mutate(diff.expressed = ifelse(padj < 0.05, "1", "0"),
                                 alpha_value = ifelse(padj < 0.05, 1, 0.01))

TMPyP4_only_plot4 = ggplot(TMPyP4_only_fc %>% arrange(diff.expressed), aes(x = log2FoldChange, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "TMPyP4-treated only TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("padj < 0.05"), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
TMPyP4_only_plot4

## grid of scatter plots
plot_grid(common_plot, nt_plot, nt_only_plot, TMPyP4_plot, TMPyP4_only_plot)

ggsave(
  glue(
    "{result_folder}TAD_bound_vs_expression_scplot.png"
  ),
  plot = last_plot(),
  width = 10,
  height = 8,
  dpi = 500,
)

plot_grid(common_plot2, nt_plot2, nt_only_plot2, TMPyP4_plot2, TMPyP4_only_plot2)

ggsave(
  glue(
    "{result_folder}TAD_bound_vs_expr-w_diff_genes.png"
  ),
  plot = last_plot(),
  width = 15,
  height = 8,
  dpi = 500,
)

ggsave(
  glue(
    "{result_folder}TAD_bound_vs_expr-w_diff_genes.pdf"
  ),
  plot = last_plot(),
  width = 15,
  height = 8
)


plot_grid(common_plot3, nt_plot3, nt_only_plot3, TMPyP4_plot3, TMPyP4_only_plot3)

ggsave(
  glue(
    "{result_folder}TAD_bound_vs_foldchanges.png"
  ),
  plot = last_plot(),
  width = 15,
  height = 8,
  dpi = 500,
)

ggsave(
  glue(
    "{result_folder}TAD_bound_vs_foldchanges.pdf"
  ),
  plot = last_plot(),
  width = 15,
  height = 8
)

plot_grid(common_plot4, nt_plot4, nt_only_plot4, TMPyP4_plot4, TMPyP4_only_plot4)

ggsave(
  glue(
    "{result_folder}TAD_bound_vs_fchanges-w_padj.png"
  ),
  plot = last_plot(),
  width = 10,
  height = 8,
  dpi = 500,
)

# examine expression levels of genes close (+/- 100 kb) to TAD boundaries
nt_boundaries = fread("../results/hi-c/NT_only_fastq_merge_200kb_annot.tsv")
nt_boundaries = nt_boundaries %>% dplyr::filter(abs(distanceToTSS) < 100000) %>%
  inner_join(., deseq2, by = c("gene_symbol" = "gene_symbol")) %>%
  dplyr::select(gene_symbol, "non-treated, rep 1" = NT_1, "non-treated, rep 2" = NT_2, 
                "TMPyP4, rep 1" = TMPyP4_1, "TMPyP4, rep 2" = TMPyP4_2,
                log2_insulation_score) %>%
  pivot_longer("non-treated, rep 1":"TMPyP4, rep 2", names_to = "condition", values_to = "norm_expr") %>% 
  mutate(hic_sample = "non-treated only") %>% 
  dplyr::select(log2_insulation_score, norm_expr, condition, hic_sample)

# create heatmaps for non-treated only and TMPyP4 only boundaries

create_heatmaps = function(path_to_boundaries,
                                     output_name) {
  
  boundaries = fread(path_to_boundaries)
  hm = boundaries %>% dplyr::filter(abs(distanceToTSS) < 100000) %>%
    inner_join(., vst, by = c("gene_symbol" = "gene_symbol")) %>%
    dplyr::select("gene_symbol", "non-treated, rep 1" = NT_1, "non-treated, rep 2" = NT_2, 
                  "TMPyP4, rep 1" = TMPyP4_1, "TMPyP4, rep 2" = TMPyP4_2)
  rows = hm %>% pull(gene_symbol)
  hm = as.matrix(hm[,2:5])
  rownames(hm) = rows
  hm = cbind(hm, rowVars(hm))
  colnames(hm)[5] <- "variance"
  hm = hm[order(hm[,5], decreasing = TRUE),][1:50,]
  hm = hm[,1:4]
  most_variable_set = rownames(hm) # top variable genes
  
  col_fun = colorRamp2(c(0, 7.5, 15), c("#440154", "#21918c", "#fde725"))
  norm_hm = Heatmap(
    hm,
    column_title = " ",
    row_title = "genes close to boundaries (+/- 100 kbp)",
    name = "norm. expr.",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    col = col_fun,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    heatmap_width = unit(4, "cm"),
    heatmap_height = unit(12, "cm"),
    show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 4),
    column_names_rot = 90
  )
  
  # fold change heatmap
  hm_fc = boundaries %>% dplyr::filter(abs(distanceToTSS) < 100000) %>%
    left_join(., fc, by = c("gene_symbol" = "gene_name")) %>%
    dplyr::filter(gene_symbol %in% most_variable_set) %>% 
    dplyr::select(gene_symbol, log2FoldChange) 
  hm_fc[is.na(hm_fc)] = 0
  rows_hm_fc = hm_fc$gene_symbol
  hm_fc = as.matrix(hm_fc[,2])
  rownames(hm_fc) = rows_hm_fc
  
  col_fun_fc = colorRamp2(c(-0.5, 0, 0.5), c("#636363", "#f0f0f0", "#feb24c"))
  hm_fc = Heatmap(
    hm_fc,
    column_title = "",
    row_title = "",
    name = "log2FC",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    col = col_fun_fc,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    heatmap_width = unit(0.5, "cm"),
    heatmap_height = unit(12, "cm"),
    show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 4),
    column_names_rot = 90,
  )
  
  # ins score heatmap
  hm_ins = boundaries %>% dplyr::filter(abs(distanceToTSS) < 100000) %>%
    left_join(., fc, by = c("gene_symbol" = "gene_name")) %>%
    dplyr::filter(gene_symbol %in% most_variable_set) %>% 
    dplyr::select(gene_symbol, log2_insulation_score) 
  hm_ins[is.na(hm_ins)] = 0
  rows_hm_ins = hm_ins$gene_symbol
  hm_ins = as.matrix(hm_ins[,2])
  rownames(hm_ins) = rows_hm_ins
  
  col_fun_ins = colorRamp2(c(0, -0.5, -1), c("#04040d", "#ad3d78", "#f8f5b6"))
  hm_ins = Heatmap(
    hm_ins,
    column_title = "",
    row_title = "",
    name = "log2InsScore",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    col = col_fun_ins,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    heatmap_width = unit(2, "cm"),
    heatmap_height = unit(12, "cm"),
    show_row_names = TRUE,
    column_names_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 4),
    column_names_rot = 90,
  )
  
  hm_list = norm_hm + hm_fc + hm_ins
  pdf(
    file = glue("{result_folder}{output_name}"),
    width = 4,
    height = 5.5
  )
  suppressWarnings(draw(hm_list, auto_adjust = FALSE))
  dev.off()
  
}

# run the heatmap function on boundaries
create_heatmaps(path_to_boundaries = "../results/hi-c/NT_only_fastq_merge_200kb_annot.tsv",
                output_name = "NT_only_heatmaps.pdf")
create_heatmaps(path_to_boundaries = "../results/hi-c/TMPyP4_only_fastq_merge_200kb_annot.tsv",
                output_name = "TMPyP4_only_heatmaps.pdf")
create_heatmaps(path_to_boundaries = "../results/hi-c/common_fastq_merge_200kb_annot.tsv",
                output_name = "common_heatmaps.pdf")


tmpyp4_boundaries = fread("../results/hi-c/TMPyP4_only_fastq_merge_200kb_annot.tsv")
tmpyp4_boundaries = tmpyp4_boundaries %>% dplyr::filter(abs(distanceToTSS) < 100000) %>%
  inner_join(., deseq2, by = c("gene_symbol" = "gene_symbol")) %>%
  dplyr::select(gene_symbol, "non-treated, rep 1" = NT_1, "non-treated, rep 2" = NT_2,
                "TMPyP4, rep 1" = TMPyP4_1, "TMPyP4, rep 2" = TMPyP4_2, log2_insulation_score) %>%
  pivot_longer("non-treated, rep 1":"TMPyP4, rep 2", names_to = "condition", values_to = "norm_expr") %>% 
  mutate(hic_sample = "TMPyP4 only") %>% 
  dplyr::select(log2_insulation_score, norm_expr, condition, hic_sample)


# violin plot of expression of genes close to TAD boundaries
vis = rbind(nt_boundaries, tmpyp4_boundaries)
  
comparisons = list(
  c("non-treated, rep 1", "TMPyP4, rep 1"),
  c("non-treated, rep 2", "TMPyP4, rep 2")
)
ggplot(vis,
       aes(x = condition, y = norm_expr, fill = hic_sample)) +
  geom_violin(trim = TRUE, color = "black") +
  scale_fill_manual(values = c("#a6bddb", "#fc9272")) +
  labs(title = "genes close to TAD boundaries (+/- 100 kb)",
       x = "RNA-Seq",
       y = "DESeq2 norm. expr.",
       fill = "TAD boundary") +
  ylim(0, 30) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(
      size = 15,
      color = "black",
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(size = 13, color = "black")
  ) + stat_compare_means(comparisons = comparisons, label = "p.signif") +
  stat_compare_means(label.y = 20, label.x = 2.5, label.sep = "\n", size = 3)


ggsave(
  glue("{result_folder}TAD_boundary_expr_violins_exp.pdf"),
  plot = last_plot(),
  width = 7,
  height = 5,
  device = "pdf"
)

ggsave(
  glue("{result_folder}TAD_boundary_expr_violins_exp.png"),
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300
)





                                   