# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("ggpubr")
})

# comparison of Daan's old mESC single-end RNA-Seq (2015) and the newly sequenced TMPyP4 RNA-Seq
# result folder
rnaseq_data = "../data/RNA-Seq/"
result_folder = "../results/rna_seq_deseq2/"

rna1 = fread(glue("{rnaseq_data}mESC_bulk_RNASeq_20151023/star_rsem/wildtype_mESC_1.genes.results"))
rna2 = fread(glue("{rnaseq_data}mESC_bulk_RNASeq_20151023/star_rsem/wildtype_mESC_2.genes.results"))
rna2015 = rna1 %>% inner_join(., rna2, by = "gene_id") %>% dplyr::select(gene_id, TPM.x, TPM.y) %>% 
  mutate(., TPM_mean = rowMeans(select(., TPM.x:TPM.y), na.rm = TRUE)) %>% dplyr::select(gene_id, TPM_mean) %>% 
  dplyr::filter(TPM_mean < 7500)

rna_tmpyp4 = fread(glue("{rnaseq_data}mESC_bulk_RNASeq_TMPyP4_study/star_rsem/rsem.merged.gene_tpm.tsv"))
rna_tmpyp4 = rna_tmpyp4 %>% mutate(., TPM_non_trt_mean = rowMeans(select(., c(NT_1, NT_2)), na.rm = TRUE)) %>% 
  mutate(., TPM_TMPyP4_mean = rowMeans(select(., c(TMPyP4_1, TMPyP4_2)), na.rm = TRUE)) %>% dplyr::select(gene_id, ends_with("mean")) %>% 
  dplyr::filter(TPM_TMPyP4_mean < 7500)

rna_all = rna_tmpyp4 %>% inner_join(., rna2015, by = c("gene_id"))

ggplot(rna_all, aes(x = TPM_non_trt_mean, y = TPM_mean)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    fill = "#fc9272"
  ) +
  # scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("end", "start")) +
  scale_color_manual(values = "black") +
  labs(
    title = "mESC RNA-Seq 2015 vs. mESC non-trt RNA-Seq",
    x = "TPM (TMPyP4 non-trt)",
    y = "TPM (mESC 2015)",
    fill = ""
  ) +
  #ylim(-5, 2) +
  guides(alpha = FALSE, size = FALSE, fill = FALSE, reverse = TRUE) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  stat_cor(method = "pearson", label.x = 2000, label.y = 6000, size = 10)

# non_tr_low = rna_all %>% dplyr::filter(TPM_non_trt_mean < 500 & TPM_mean > 2000)

ggsave(
  plot = last_plot(),
  glue("{result_folder}mESC_RNA-Seqs_scatter.pdf"),
  width = 8,
  height = 6
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}mESC_RNA-Seqs_scatter.png"),
  width = 8,
  height = 6,
  dpi = 300
)


