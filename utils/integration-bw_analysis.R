suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("glue")
  library("ggpubr")
  library("ggrepel")
  library("ggrastr")
})

# result folder
result_folder = "../results/rna_seq_deseq2/"

# RNA-Seq data
deseq2 = "../results/rna_seq_deseq2/reg_log_norm_counts.tsv"
vst = "../results/rna_seq_deseq2/vst_norm_counts.tsv"
deseq2 = fread(deseq2)
vst = fread(vst)

fc = fread("../results/rna_seq_deseq2/RNA-Seq_treat_vs_contr_DESeq2_fc.tsv")
sign = fc %>% filter(abs(log2FoldChange) > 0.25) %>% pull(gene_name)

# let's see the TMPyP4 vs. control RNA-Seq fold changes of those genes that are close to TSS where CTCF enriched upon AID (6 hours) treatment
# volcano plot by ggplot2
create_volc = function(bw_comp_output1,
                       bw_comp_output2,
                       fc_filter1 = 2,
                       fc_filter2 = 5,
                       distance_filter = 3000,
                       label1,
                       label2,
                       volcano_title,
                       output_file,
                       export = TRUE) {
  
  # TMPyP4 treatment
  # bw_loci - G4 signals following 6h treatment over CTCF WT peaks
  bw_comp_output1 = fread(bw_comp_output1)
  # bw_loci - CTCFc signals following 6h treatment over CTCF WT peaks
  bw_comp_output2 = fread(bw_comp_output2)
  
  joined = inner_join(
    bw_comp_output1,
    bw_comp_output2,
    by = "gene_symbol",
    suffix = c("_1", "_2")
  )
  
  joined_filt = joined %>% dplyr::filter(
    fold_change_1 > fc_filter1 &
      fold_change_2 > fc_filter2 &
      abs(distanceToTSS_1) < distance_filter & abs(distanceToTSS_2)
  ) %>%
    left_join(., fc, by = c("gene_symbol" = "gene_name"))

  volc = joined_filt %>% drop_na() %>% mutate(
    group = case_when(
      log2FoldChange > 0.10 & pvalue < 0.05 ~ "up",
      log2FoldChange < -0.10 & pvalue < 0.05 ~ "down",
      log2FoldChange >= -0.10 & log2FoldChange <= 0.10 ~ "unaltered"
    )
  ) %>%
    mutate(
      sign_label = case_when(
        log2FoldChange > 0.10 & -log10(pvalue) > 3 ~ gene_symbol,
        log2FoldChange < -0.10 &
          -log10(pvalue) > 3 ~ gene_symbol,
        log2FoldChange >= -0.10 &
          log2FoldChange <= 0.10 ~ ""
      )
    )
  labels = volc %>% pull(sign_label)
  
  cols = c("up" = "#fc9272",
           "down" = "#a1d99b",
           "unaltered" = "grey")
  sizes = c("up" = 4,
            "down" = 4,
            "unaltered" = 2)
  alphas = c("up" = 1,
             "down" = 1,
             "unaltered" = 0.5)
  
  ggplot_volc = volc %>%
    ggplot(aes(
      x = log2FoldChange,
      y = -log10(pvalue),
      fill = group,
      size = group,
      alpha = group
    )) +
    ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters
                             colour = "black") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") +
    geom_vline(xintercept = c(-0.10, 0.10),
               linetype = "dashed") +
    scale_fill_manual(values = cols) + # Modify point colour
    scale_size_manual(values = sizes) + # Modify point size
    scale_alpha_manual(values = alphas) + # Modify point transparency
    scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),
                       limits = c(-1.5, 1.5)) +
    labs(title = glue("TMPyP4 vs. ctrl RNA-Seq, regions: {volcano_title}"),
         subtitle = glue("cut offs: |log2FC| > 0.10, p-value < 0.05, {label1} FC > {fc_filter1}, {label2} FC > {fc_filter2}, TSS +/- {distance_filter/1000} kbp"),
         x = "log2FoldChange",
         y = "-log10 p-value",
         fill = " ") +
    guides(alpha = "none",
           size = "none",
           fill = guide_legend(override.aes = list(size = 5))) +
    theme_minimal() +
    theme(
      text = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 18, color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.title.x = element_text(size = 23, color = "black"),
      axis.title.y = element_text(size = 23, color = "black")
    ) +
    geom_text_repel(label = labels, size = 6, max.overlaps = 20)
  print(ggplot_volc)
  
  if (export) {
    ggsave(
      glue("{result_folder}{output_file}.pdf"),
      plot = ggplot_volc,
      width = 10,
      height = 10,
      device = "pdf"
    )
    
    ggsave(
      glue("{result_folder}{output_file}.png"),
      plot = ggplot_volc,
      width = 10,
      height = 10,
      dpi = 300
    )
    
  }
  
  return(ggplot_volc)
  
}

create_volc(bw_comp_output1 = "../results/cutntag/G4down_AID_6h_vs_0h_at_CTCF_peaks.tsv",
            bw_comp_output2 = "../results/cutntag/CTCFup_AID_6h_vs_0h_at_CTCF_peaks.tsv",
            volcano_title = "G4 lost, CTCF enriched upon AID (6h)",
            label1 = "G4 down",
            label2 = "CTCF up",
            output_file = "G4down_AID_6h_CTCFup_AID_6h_RNA_Seq_volc")

create_volc(bw_comp_output1 = "../results/cutntag/G4up_AID_6h_vs_0h_at_CTCF_peaks.tsv",
            bw_comp_output2 = "../results/cutntag/CTCFdown_AID_6h_vs_0h_at_CTCF_peaks.tsv",
            volcano_title = "G4 gained, CTCF lost upon AID (6h)",
            label1 = "G4 up",
            label2 = "CTCF down",
            output_file = "G4up_AID_6h_CTCFdown_AID_6h_RNA_Seq_volc")
                





