suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("ComplexHeatmap")
  library("circlize")
  library("cowplot")
})

result_folder = "../results/deeptools/"

correlations = fread("../results/cutntag/Wulfridge_CTCF_ChIPSeq_G4_no_CTCF_proms_correlations.tab")

correlations$V1

# samples = c("NPC G4 bulk C&T", "MEF G4 bulk C&T", "mESC G4 bulk C&T", 
#             "MEF G4 scC&T (cluster 0)",  "mESC G4 scC&T (cluster 1)", 
#             "brain cells G4 scC&T (cluster 1)",
#             "oligodendr. G4 scC&T (cluster 0)")

samples = c("Wulfridge CTCF ChIP-Seq, PDS",
            "Wulfridge CTCF ChIP-Seq, mock",
            "Wulfridge CTCF ChIP-Seq, PhenDC3",
            "G4 Cut&Tag, non-trt.", "G4 Cut&Tag, TMPyP4 2h",
            "G4 Cut&Tag, TMPyP4 6h", "CTCF Cut&Tag, TMPyP4 0h",
            "CTCF Cut&Tag, TMPyP4 6h")


correlations = correlations %>% dplyr::select(-V1)
correlations = as.matrix(correlations)
rownames(correlations) = samples
colnames(correlations) = samples


col_fun = colorRamp2(c(-1, 0, 1), c("#9ecae1", "white", "#fc9272"))

# png(
#   file = glue("{result_folder}Pearson_bw_corr_hm.png"),
#   width = 8,
#   height = 8,
#   units = 'in',
#   res = 500
# )
pdf(
  file = glue("{result_folder}Wulfridge_CTCF_ChIP-Seq-G4only_proms_bw_corr_hm.pdf"),
  width = 8,
  height = 8
)
corr_hm = Heatmap(
  correlations,
  name = "Pearson",
  clustering_distance_rows = "pearson",
  col = col_fun,
  row_title = "",
  column_title = "G4 only promoters",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", correlations[i, j]), x, y, gp = gpar(fontsize = 10))
  },
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  width = unit(7, "cm"),
  height = unit(7, "cm"),
)
corr_hm

dev.off()


png(
  file = glue("{result_folder}Wulfridge_CTCF_ChIP-Seq-G4only_proms_bw_corr_hm.png"),
  width = 16,
  height = 16,
  unit = "cm",
  res = 300
)
corr_hm = Heatmap(
  correlations,
  name = "Pearson",
  clustering_distance_rows = "pearson",
  col = col_fun,
  row_title = "",
  column_title = "G4 only promoters",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", correlations[i, j]), x, y, gp = gpar(fontsize = 10))
  },
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  width = unit(7, "cm"),
  height = unit(7, "cm"),
)
corr_hm
dev.off()



