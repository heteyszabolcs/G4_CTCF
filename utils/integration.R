# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("cowplot")
})

result_folder = "../results/hi-c/"

tad_annot = "../results/hi-c/"
tad_annot = list.files(tad_annot, pattern = 'TAD_annots.tsv', full.names = TRUE)

tad = "../data/Hi-C/bed/"
tad = list.files(tad, pattern = "bed", full.names = TRUE)

deseq2 = "../results/rna_seq_deseq2/reg_log_norm_counts.tsv"
deseq2 = fread(deseq2)

fc = fread("../results/rna_seq_deseq2/RNA-Seq_treat_vs_contr_DESeq2_fc.tsv")
sign = fc %>% filter(abs(log2FoldChange) > 0.5) %>% pull(gene_name)

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

# highlight genes with higher than abs(0.5) DESeq2 fold change
common_long = common_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "1", "0"),
                                     alpha_value = ifelse(gene_symbol %in% sign, 1, 0.01))

common_plot2 = ggplot(common_long %>% arrange(diff.expressed), aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "common TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("diff. expr."), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
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
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "common TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  xlim(-1.5, 1.5) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
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

nt_long = nt_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "1", "0"),
                                     alpha_value = ifelse(gene_symbol %in% sign, 1, 0.01))

nt_plot2 = ggplot(nt_long %>% arrange(diff.expressed), aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "non-treated TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("diff. expr."), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
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
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "non-treated TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  xlim(-1.5, 1.5) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
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

nt_only_long = nt_only_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "1", "0"),
                             alpha_value = ifelse(gene_symbol %in% sign, 1, 0.01))

nt_only_plot2 = ggplot(nt_only_long %>% arrange(diff.expressed), aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "non-treated only TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("diff. expr."), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
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
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "non-treated TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  xlim(-1.5, 1.5) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
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

TMPyP4_long = TMPyP4_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "1", "0"),
                                       alpha_value = ifelse(gene_symbol %in% sign, 1, 0.01))

TMPyP4_plot2 = ggplot(TMPyP4_long %>% arrange(diff.expressed), aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "TMPyP4-treated TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("diff. expr."), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
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
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "TMPyP4-treated TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  xlim(-1.5, 1.5) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
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

TMPyP4_only_long = TMPyP4_only_long %>% mutate(diff.expressed = ifelse(gene_symbol %in% sign, "1", "0"),
                                     alpha_value = ifelse(gene_symbol %in% sign, 1, 0.01))

TMPyP4_only_plot2 = ggplot(TMPyP4_only_long %>% arrange(diff.expressed), aes(x = NT_1, y = log2_ins_score)) +
  geom_point(size = 2, aes(colour = diff.expressed, alpha = alpha_value)) +
  scale_color_manual(values = c("#f0f0f0", "#f03b20")) +
  labs(
    title = "TMPyP4-treated only TAD boundaries",
    x = "DESeq2 norm. expr.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  guides(colour = guide_legend("diff. expr."), alpha = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
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
  geom_point(size = 2, aes(colour = ins_score_type)) +
  # geom_point(fill = "black", size = 7, shape = 1) +
  scale_color_manual(values = c("#fdbb84", "#a6bddb")) +
  labs(
    title = "TMPyP4-treated only TAD boundaries",
    x = "DESeq2 log2 fold change.",
    y = "log2 insulation score",
    color = " "
  ) +
  ylim(-5, 2) +
  xlim(-1.5, 1.5) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )
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

## sum up
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
  width = 10,
  height = 8,
  dpi = 500,
)

plot_grid(common_plot3, nt_plot3, nt_only_plot3, TMPyP4_plot3, TMPyP4_only_plot3)

ggsave(
  glue(
    "{result_folder}TAD_bound_vs_foldchanges.png"
  ),
  plot = last_plot(),
  width = 10,
  height = 8,
  dpi = 500,
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
