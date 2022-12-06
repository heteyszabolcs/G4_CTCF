# packages
suppressPackageStartupMessages({
  library("wigglescout")
  library("data.table")
  library("tidyverse")
  library("glue")
  library("ggplot2")
  library("glue")
  library("ggfortify")
  library("ggrepel")
  library("ggpubr")
})

# folder
bigwig_folder = "../data/CutNTag/bw/"
rnaseq_folder = "../results/rna_seq_deseq2/"
result_folder = "../results/cutntag/"

# Cut&Tag bigwigs
bigwigs = list.files(bigwig_folder, pattern = "*.norm.bw", full.names = TRUE)

# gene regions (+ 500 bp from TSS and gene body)
gene_regions = list.files(rnaseq_folder, pattern = "*_500bp_upstream.bed", full.names = TRUE)
high = "../results/rna_seq_deseq2/NT_high_500bp_upstream.bed"
low = "../results/rna_seq_deseq2/NT_low_500bp_upstream.bed"
lowmedium = "../results/rna_seq_deseq2/NT_lowmedium_500bp_upstream.bed"
highmedium = "../results/rna_seq_deseq2/NT_highmedium_500bp_upstream.bed"

# PCA function
plot_pca = function(path_to_bed = lowmedium, title) {
  aggr = bw_loci(bigwigs, loci = path_to_bed)
  aggr = as.data.frame(aggr)
  aggr = aggr %>% dplyr::select(contains("norm")) %>% drop_na() %>% dplyr::select(
    "CTCF 0h AID" = "CTCF_AID_0h_merge_5million.norm",
    "CTCF 6h AID" = "CTCF_AID_6h_merge.norm",
    "CTCF TMPyP4 0h" = "CTCF_TMPyP4_0h_merge_5million.norm",
    "CTCF TMPyP4 6h" = "CTCF_TMPyP4_6h_merge_5million.norm",
    "G4 0h AID" = "G4_0h_AID_5million.norm",
    "G4 24h AID" = "G4_24h_AID_5million.norm",
    "G4 48h AID" = "G4_48h_AID_5million.norm",
    "G4 6h AID" = "G4_6h_AID_5million.norm"
  )
  aggr = as.matrix(aggr)
  
  # generate PCA
  pca = prcomp(t(aggr), scale. = FALSE)
  # add colors
  values = c(
    "#a6bddb",
    "#2b8cbe",
    "#a1d99b",
    "#31a354",
    "#fee0d2",
    "#fc9272",
    "#de2d26",
    "#521815"
  )
  # plot PCA
  p = autoplot(pca,
               data = t(aggr),
               loadings = FALSE,
               label = FALSE) +
    geom_point(
      size = 12,
      shape = 21,
      colour = "black",
      aes(fill = factor(
        rownames(pca$x),
        levels = c(
          "CTCF 0h AID",
          "CTCF 6h AID",
          "CTCF TMPyP4 0h",
          "CTCF TMPyP4 6h",
          "G4 0h AID",
          "G4 24h AID",
          "G4 48h AID",
          "G4 6h AID"
        )
      ))
    ) +
    labs(title = title,
         subtitle = "window: 500 bp upsteam to TSS - TES",
         color = "",
         fill = "") +
    ylim(-1, 1) +
    xlim(-1, 1) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 8),
      axis.text.x = element_text(size = 14, color = "black")
    ) +
    scale_fill_manual(values = values)
  
  return(print(p))
  
}

# generate PCA plots
high_pca = plot_pca(high, title = "PCA on G4/CTCF signals of high expressed genes")
mediumhigh_plot = plot_pca(highmedium, title = "PCA on G4/CTCF signals of medium-high expressed genes")
mediumlow_pca = plot_pca(lowmedium, title = "PCA on G4/CTCF signals of medium-low expressed genes")
low_pca = plot_pca(low, title = "PCA on G4/CTCF signals of low expressed genes")

pca_list = list(high_pca, mediumhigh_plot, mediumlow_pca, low_pca)
ggarrange(plotlist = pca_list)

# export
ggsave(
  plot = last_plot(),
  glue("{result_folder}pca_plots_on_gene_groups.pdf"),
  width = 12,
  height = 10
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}pca_plots_on_gene_groups.png"),
  width = 12,
  height = 10,
  dpi = 300
)
