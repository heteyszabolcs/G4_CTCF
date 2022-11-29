# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("ComplexHeatmap")
  library("circlize")
  library("ggrepel")
  library("ggpubr")
  library("glue")
  library("caret")
  library("caTools")
})

set.seed(42)

# result folder
result_folder = "../results/hi-c/"
cnt_folder = "../results/cutntag/"

# RNA-Seq data
deseq2 = "../results/rna_seq_deseq2/reg_log_norm_counts.tsv"
vst = "../results/rna_seq_deseq2/vst_norm_counts.tsv"
deseq2 = fread(deseq2)
vst = fread(vst)

fc = fread("../results/rna_seq_deseq2/RNA-Seq_treat_vs_contr_DESeq2_fc.tsv")
sign = fc %>% filter(abs(log2FoldChange) > 0.25) %>% pull(gene_name)


create_ins_fc_summary = function(annot_bound_path, k = 3, type = "common") {
  
  boundary_data = fread(annot_bound_path)
  km = kmeans(boundary_data$log2_insulation_score, k, nstart = 25)
  boundary_data = boundary_data %>% mutate(kmeans = as.character(km$cluster))
  
  # boxplot - kmeans 
  km_plot = ggplot(boundary_data, aes(x = kmeans, y = log2_insulation_score, fill = kmeans)) +
    geom_boxplot(color = "black") +
    scale_fill_manual(values = c("#fc9272", "#9ecae1", "#a1d99b")) +
    ylim(-4, 1) +
    labs(
      title = glue("kmeans - {type} boundaries"),
      x = "",
      y = "log2InsScore",
      fill = ""
    ) +
    guides(fill = "none") +
    theme_bw() +
    theme(
      text = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 20, color = "black"),
      axis.title.y = element_text(size = 20, color = "black")
    ) +                                      
    annotate("text",
             x = 1:length(table(boundary_data$kmeans)),
             y = aggregate(log2_insulation_score ~ kmeans, boundary_data, max)[ , 2],
             label = paste("n =", as.character(table(boundary_data$kmeans)), sep = " "),
             col = "black",
             size = 7,
             vjust = - 1)
  
  # add expr features + labels
  boundary_data = boundary_data %>% 
    inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
    mutate(label = case_when(log2_insulation_score < -2 ~ gene_symbol,
                             abs(log2FoldChange) > 0.5 ~ gene_symbol)) %>% 
    mutate(close = case_when(abs(distanceToTSS) < 100000 ~ "+/- 100 kb", TRUE ~ "distant")) %>% 
    mutate(sign = case_when(gene_symbol %in% sign ~ "|log2FC| > 0.25", TRUE ~ "non-altered")) %>% 
    mutate(sign_label = case_when(sign == "|log2FC| > 0.25" ~ gene_symbol, TRUE ~ ""))
  
  # scatterplot 2 - kmeans labeled with remarkable scores 
  labeled_plot = ggplot(boundary_data, aes(x = log2FoldChange, y = log2_insulation_score)) +
    ggrastr::geom_point_rast(size = 2,
                             shape = 21,
                             colour = "black",
                             aes(fill = kmeans)) +
    scale_fill_manual(values = c("#fc9272", "#9ecae1", "#a1d99b")) +
    scale_color_manual(values = "black") +
    labs(title = glue("{type} boundaries"),
         x = "log2FoldChange",
         y = "log2InsScore",
         fill = "cluster") +
    scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),
                       limits = c(-1.5, 1.5)) +
    ylim(-4, 1) +
    guides(
      alpha = "none",
      size = "none",
      fill = guide_legend(override.aes = list(size = 5), reverse = FALSE)
    ) +
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
    geom_text_repel(aes(label = label), size = 5, max.overlaps = 50) 
  
  # scatterplot 3 - label genes close to boundaries
  distance_plot = ggplot(boundary_data, aes(x = log2FoldChange, y = log2_insulation_score)) +
    ggrastr::geom_point_rast(size = 2,
                             shape = 21,
                             colour = "black",
                             aes(fill = close)) +
    scale_fill_manual(values = c("#f03b20", "#f0f0f0")) +
    scale_color_manual(values = "black") +
    labs(title = glue("{type} boundaries"),
         x = "log2FoldChange",
         y = "log2InsScore",
         fill = "Distance to TAD") +
    scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),
                       limits = c(-1.5, 1.5)) +
    ylim(-4, 1) +
    guides(
      alpha = "none",
      size = "none",
      fill = guide_legend(override.aes = list(size = 5), reverse = FALSE)
    ) +
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
    geom_text_repel(aes(label = label), size = 5, max.overlaps = 80)
  
  # scatterplot 4 - labeled of altered genes
  fc_plot = ggplot(boundary_data, aes(x = log2FoldChange, y = log2_insulation_score)) +
    ggrastr::geom_point_rast(size = 2,
                             shape = 21,
                             colour = "black",
                             aes(fill = sign)) +
    scale_fill_manual(values = c("#2c7fb8", "#f0f0f0")) +
    scale_color_manual(values = "black") +
    labs(title = glue("{type} boundaries"),
         x = "log2FoldChange",
         y = "log2InsScore",
         fill = "Expr. status") +
    scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),
                       limits = c(-1.5, 1.5)) +
    ylim(-4, 1) +
    guides(
      alpha = "none",
      size = "none",
      fill = guide_legend(override.aes = list(size = 5), reverse = FALSE)
    ) +
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
    geom_text_repel(aes(label = sign_label), size = 5, max.overlaps = 50)

  # arrange
  panel = ggarrange(km_plot, labeled_plot, fc_plot, distance_plot, ncol = 2, nrow = 2)

  # export
  ggsave(
    glue("{result_folder}{type}_boundary_km_scatter_sum.png"),
    plot = panel,
    width = 12,
    height = 10,
    dpi = 500,
  )
  ggsave(
    glue("{result_folder}{type}_boundary_km_scatter_sum.pdf"),
    plot = panel,
    width = 12,
    height = 10,
    device = "pdf"
  )
  
  return(print(panel))
  
}

create_ins_fc_summary(annot_bound_path = "../data/Hi-C/bed/common_fastq_merge_200kb_annot.tsv", type = "common")
create_ins_fc_summary(annot_bound_path = "../data/Hi-C/bed/NT_only_fastq_merge_200kb_annot.tsv", type = "non-trt only")
create_ins_fc_summary(annot_bound_path = "../data/Hi-C/bed/TMPyP4_only_fastq_merge_200kb_annot.tsv", type = "TMPyP4 only")


# creating statistical model to examine if we can predict insulation scores by expressional changes
# insulation scores are subset in terms of conditions
# package: caret
common = fread("../data/Hi-C/bed/common_fastq_merge_200kb_annot.tsv") %>% 
  mutate(type = "common") %>% 
  mutate(region = paste(seqnames, start, end, sep = "_"))
nt_only = fread("../data/Hi-C/bed/NT_only_fastq_merge_200kb_annot.tsv") %>% 
  mutate(type = "non-treated only") %>% 
  mutate(region = paste(seqnames, start, end, sep = "_"))
tmpyp4_only = fread("../data/Hi-C/bed/TMPyP4_only_fastq_merge_200kb_annot.tsv") %>% 
  mutate(type = "TMPyP4 only") %>% 
  mutate(region = paste(seqnames, start, end, sep = "_"))

# train test split function (source: train_test_split.R)
train_test_split <- function(df) {
  set.seed(42)
  df = as.data.frame(df)
  sample = sample.split(df, SplitRatio = 0.8)
  train = subset(df, sample == TRUE)
  test  = subset(df, sample == FALSE)
  return(list(train, test))
}

# linear regression on common boundaries
common_with_fc = common %>% inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  dplyr::filter(abs(distanceToTSS) < 10000) %>% 
  mutate(strong_insulator = ifelse(log2_insulation_score < -0.5, 1, 0)) %>% 
  mutate(strong_insulator = as.factor(strong_insulator))

common_with_fc_train <- train_test_split(common_with_fc)[[1]]
common_with_fc_test <- train_test_split(common_with_fc)[[2]]

lm_common = train(log2FoldChange ~ log2_insulation_score,
                   data = common_with_fc_train,
                   method = "lm")

common_rmse = round(lm_common$results$RMSE, 2)
common_rsquared = round(lm_common$results$Rsquared, 2)
common_mae = round(lm_common$results$MAE, 2)

lm_common_plot = ggplot(common_with_fc, aes(y = log2FoldChange, x = log2_insulation_score)) + 
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = "type")
  ) +
  scale_fill_manual(values = "red") +
  labs(
    title = "common lm: insulation score ~´expression of genes +/- 100 kb TSS",
    x = "log2InsScore",
    y = "log2FoldChange",
    fill = ""
  ) +
  scale_x_continuous(breaks = c(seq(-1.5, 0.5, 0.5)),       
                     limits = c(-1.5, 0.5)) +
  ylim(-5, 5) +
  ylim(-1, 1) +
  geom_smooth(method = "lm", colour = "black") +
  guides(alpha = "none", size = "none", fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) + 
  geom_text(aes(label = c(glue("RMSE: {common_rmse}, R squared: {common_rsquared}, MAE: {common_mae}"))), size = 4, x = -0.5, y = 0.70)
lm_common_plot


# linear regression on non-treated only boundaries
nt_only_with_fc = nt_only %>% inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  dplyr::filter(abs(distanceToTSS) < 10000) %>% 
  mutate(strong_insulator = ifelse(log2_insulation_score < -0.5, 1, 0)) %>% 
  mutate(strong_insulator = as.factor(strong_insulator))

nt_only_with_fc_train <- train_test_split(nt_only_with_fc)[[1]]
nt_only_with_fc_test <- train_test_split(nt_only_with_fc)[[2]]

nt_only_glm_model = train(strong_insulator ~ log2FoldChange, 
                         data = nt_only_with_fc_train, 
                         method = "glm",
                         family = "binomial") # binomial defines the logistic regression
summary(nt_only_glm_model)

lm_nt_only = train(log2FoldChange ~ log2_insulation_score,
                  data = nt_only_with_fc_train,
                  method = "lm")

nt_only_rmse = round(lm_nt_only$results$RMSE, 2)
nt_only_rsquared = round(lm_nt_only$results$Rsquared, 2)
nt_only_mae = round(lm_nt_only$results$MAE, 2)

lm_nt_only_plot = ggplot(nt_only_with_fc, aes(y = log2FoldChange, x = log2_insulation_score)) + 
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = "type")
  ) +
  scale_fill_manual(values = "red") +
  labs(
    title = "non-trt only lm: insulation score ~´expression of genes +/- 100 kb TSS",
    x = "log2InsScore",
    y = "log2FoldChange",
    fill = ""
  ) +
  scale_x_continuous(breaks = c(seq(-1.5, 0.5, 0.5)),       
                     limits = c(-1.5, 0.5)) +
  ylim(-1, 1) +
  xlim(-1.5, 1) +
  geom_smooth(method = "lm", colour = "black") +
  guides(alpha = "none", size = "none", fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )  + 
  geom_text(aes(label = c(glue("RMSE: {nt_only_rmse}, R squared: {nt_only_rsquared}, MAE: {nt_only_mae}"))), size = 4, x = -0.5, y = 0.70)
lm_nt_only_plot

# linear regression on TMPyP4 only boundaries
tmpyp4_only_with_fc = tmpyp4_only %>% inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  dplyr::filter(abs(distanceToTSS) < 10000) %>% 
  mutate(strong_insulator = ifelse(log2_insulation_score < -0.5, 1, 0)) %>% 
  mutate(strong_insulator = as.factor(strong_insulator))

tmpyp4_only_with_fc_train <- train_test_split(tmpyp4_only_with_fc)[[1]]
tmpyp4_only_with_fc_test <- train_test_split(tmpyp4_only_with_fc)[[2]]

tmpyp4_only_glm_model = train(strong_insulator ~ log2FoldChange, 
                          data = tmpyp4_only_with_fc_train, 
                          method = "glm",
                          family = "binomial") # binomial defines the logistic regression
summary(tmpyp4_only_glm_model)

lm_tmpyp4_only = train(log2FoldChange ~ log2_insulation_score,
                   data = tmpyp4_only_with_fc_train,
                   method = "lm")

tmpyp4_only_rmse = round(lm_tmpyp4_only$results$RMSE, 2)
tmpyp4_only_rsquared = round(lm_tmpyp4_only$results$Rsquared, 2)
tmpyp4_only_mae = round(lm_tmpyp4_only$results$MAE, 2)
  
tmpyp4_only_plot = ggplot(tmpyp4_only_with_fc, aes(y = log2FoldChange, x = log2_insulation_score)) + 
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = "type")
  ) +
  scale_fill_manual(values = "red") +
  labs(
    title = "TMPyP4 only lm: insulation score ~´expression of genes +/- 100 kb TSS",
    x = "log2InsScore",
    y = "log2FoldChange",
    fill = ""
  ) +
  scale_x_continuous(breaks = c(seq(-1.5, 0.5, 0.5)),       
                     limits = c(-1.5, 0.5)) +
  ylim(-1, 1) +
  xlim(-1.5, 1) +
  geom_smooth(method = "lm", colour = "black") +
  guides(alpha = "none", size = "none", fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) + 
  geom_text(aes(label = c(glue("RMSE: {tmpyp4_only_rmse}, R squared: {tmpyp4_only_rsquared}, MAE: {tmpyp4_only_mae}"))), size = 4, x = -0.5, y = 0.70)
tmpyp4_only_plot

lm_plots = ggarrange(lm_nt_only_plot, lm_common_plot, tmpyp4_only_plot)

# export
ggsave(
  glue("{result_folder}insulation_sc_vs_foldchange-lm_plots.png"),
  plot = lm_plots,
  width = 12,
  height = 10,
  dpi = 500,
)
ggsave(
  glue("{result_folder}insulation_sc_vs_foldchange-lm_plots.pdf"),
  plot = lm_plots,
  width = 12,
  height = 10,
  device = "pdf"
)

# creating statistical model to examine if log2FoldChanges can be predicted by sites with elevated G4 level and lost CTCF level. 
# fold changes of G4 and CTCF come from bw_loci analysis of bw signals (source: bw_comparison)
# package: caret
g4gained_peaks = fread("../results/cutntag/G4up_AID_6h_vs_0h_at_CTCF_peaks.tsv")
ctcflost_peaks = fread("../results/cutntag/CTCFdown_AID_6h_vs_0h_at_CTCF_peaks.tsv")

joined = inner_join(
  g4gained_peaks,
  ctcflost_peaks,
  by = "gene_symbol",
  suffix = c("_1", "_2")
)

joined = joined %>% dplyr::filter(
  fold_change_1 > 2 & # G4 up
    fold_change_1 < 30 &
    fold_change_2 > 5 & # CTCF down
    abs(distanceToTSS_1) < 3000 & abs(distanceToTSS_2) < 3000
) %>%
  left_join(., fc, by = c("gene_symbol" = "gene_name"))

joined = joined %>% drop_na() %>% mutate(
  group = case_when(
    log2FoldChange > 0.10 & pvalue < 0.05 ~ "up", 
    log2FoldChange < -0.10 & pvalue < 0.05 ~ "down",
    log2FoldChange >= -0.10 & log2FoldChange <= 0.10 ~ "unaltered"
  )
) %>% 
  dplyr::select(G4_fc = fold_change_1, CTCF_fc = fold_change_2, log2FoldChange, gene_symbol)

joined = drop_na(joined)

joined_train <- train_test_split(joined)[[1]]
joined_test <- train_test_split(joined)[[2]]

# linear model: expression change ~ enriched G4 signal (and lost CTCF signal)
lm_g4gained_fc = train(log2FoldChange ~ G4_fc,
                               data = joined_train,
                               method = "lm")
summary(lm_g4gained_fc)

g4gained_rmse = round(lm_g4gained_fc$results$RMSE, 2)
g4gained_rsquared = round(lm_g4gained_fc$results$Rsquared, 2)
g4gained_mae = round(lm_g4gained_fc$results$MAE, 2)

g4gained_fc_plot = ggplot(joined_train, aes(y = log2FoldChange, x = G4_fc)) + 
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = "type")
  ) +
  scale_fill_manual(values = "red") +
  labs(
    title = "Gained G4 regions: log2FC ~ peaks gained G4, lost CTCF (+/- 3 kb to TSS)",
    x = "G4 fold change (6h AID to 0h AID)",
    y = "RNA-Seq log2FoldChange",
    fill = ""
  ) +
  # scale_x_continuous(breaks = c(seq(-1.5, 0.5, 0.5)),       
  #                    limits = c(-1.5, 0.5)) +
  # ylim(-1, 1) +
  # xlim(-1.5, 1) +
  geom_smooth(method = "lm", colour = "black") +
  guides(alpha = "none", size = "none", fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) + 
  geom_text(aes(label = c(glue("RMSE: {g4gained_rmse}, R squared: {g4gained_rsquared}, MAE: {g4gained_mae}"))), size = 4, x = 20, y = 0.23)
g4gained_fc_plot

# export
ggsave(
  glue("{cnt_folder}G4enriched_vs_RNASeqfoldchange-lm_plot.png"),
  plot = g4gained_fc_plot,
  width = 8,
  height = 8,
  dpi = 500,
)
ggsave(
  glue("{cnt_folder}G4enriched_vs_RNASeqfoldchange-lm_plot.pdf"),
  plot = g4gained_fc_plot,
  width = 8,
  height = 8,
  device = "pdf"
)
