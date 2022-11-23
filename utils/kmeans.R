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
})

set.seed(42)

# result folder
result_folder = "../results/hi-c/"

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











