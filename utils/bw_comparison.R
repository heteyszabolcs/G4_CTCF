suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("wigglescout")
  library("GenomicRanges")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
})

# annotation script
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# result folder
result_folder = "../results/cutntag/"

# folders
bw = "../data/CutNTag/bw/"
peak = "../data/CutNTag/bed/"


calc_bigwig_fc = function(bigwig,
                          bg_bigwig,
                          subset,
                          output_file) {
  # wigglescout's bw_loci
  fc = bw_loci(bigwig, bg_bwfiles = bg_bigwig, loci = subset)
  fc = as.data.frame(fc)
  fc = fc[!is.infinite(fc[, 6]),]
  fc = fc[!is.nan(fc[, 6]), ]
  
  feature = colnames(fc)[6]
  
  # annotation
  fc = mm10_annotation(
    regions = fc,
    start_col = "start",
    end_col = "end",
    seqname_col = "seqnames",
    feature_1 = feature,
    feature_2 = feature
  )
  
  fc = fc %>%
    dplyr::select(
      seqnames,
      start,
      end,
      fold_change = feature_1,
      gene_symbol = SYMBOL,
      distanceToTSS
    ) %>%
    dplyr::filter(fold_change >= 2)
  
  write_tsv(fc, glue("{result_folder}{output_file}"))
  
  return(fc)
}

calc_bigwig_fc(bigwig = glue("{bw}CTCF_AID_0h_merge_5million.norm.bw"), 
               bg_bigwig = glue("{bw}G4_0h_AID_5million.norm.bw"),
               subset = glue("{peak}CTCF_WT_peaks.bed"),
               output_file = "CTCF_vs_G4_0H_AID.tsv")

# AID treatment
ctcf_aid_ctcf_increasing = calc_bigwig_fc(bigwig = glue("{bw}CTCF_AID_6h_merge.norm.bw"), 
               bg_bigwig = glue("{bw}CTCF_AID_0h_merge_5million.norm.bw"),
               subset = glue("{peak}CTCF_WT_peaks.bed"),
               output_file = "CTCFup_AID_6h_vs_0h_at_CTCF_peaks.tsv")

ctcf_aid_ctcf_decreasing = calc_bigwig_fc(bigwig = glue("{bw}CTCF_AID_0h_merge_5million.norm.bw"), 
                          bg_bigwig = glue("{bw}CTCF_AID_6h_merge.norm.bw"),
                          subset = glue("{peak}CTCF_WT_peaks.bed"),
                          output_file = "CTCFdown_AID_6h_vs_0h_at_CTCF_peaks.tsv")


g4_aid_g4_increasing = calc_bigwig_fc(bigwig = glue("{bw}G4_6h_AID_5million.norm.bw"), 
                             bg_bigwig = glue("{bw}G4_0h_AID_5million.norm.bw"),
                             subset = glue("{peak}CTCF_WT_peaks.bed"),
                             output_file = "G4up_AID_6h_vs_0h_at_CTCF_peaks.tsv")

g4_aid_g4_decreasing = calc_bigwig_fc(bigwig = glue("{bw}G4_0h_AID_5million.norm.bw"), 
                                      bg_bigwig = glue("{bw}G4_6h_AID_5million.norm.bw"),
                                      subset = glue("{peak}CTCF_WT_peaks.bed"),
                                      output_file = "G4down_AID_6h_vs_0h_at_CTCF_peaks.tsv")




