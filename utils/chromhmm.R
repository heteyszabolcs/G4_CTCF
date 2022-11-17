suppressPackageStartupMessages({
  library("tidyverse")
  library("wigglescout")
  library("data.table")
  library("ggplot2")
  library("glue")
})

# folders
result_folder = "../results/cutntag/"
bigwig_folder = "../data/CutNTag/bw/"

bigwigs = list.files(bigwig_folder, full.names = TRUE)
chromhmm = plot_bw_loci_summary_heatmap(bigwigs, loci = "../data/ESC_ChromHMM15_mm10.bed")

ggsave(
  glue("{result_folder}wigglescout_ChromHMM15_mm10_analysis.pdf"),
  plot = chromhmm,
  width = 10,
  height = 10,
  device = "pdf"
)

ggsave(
  glue("{result_folder}wigglescout_ChromHMM15_mm10_analysis.png"),
  plot = chromhmm,
  width = 10,
  height = 10,
  rbi = 300
)