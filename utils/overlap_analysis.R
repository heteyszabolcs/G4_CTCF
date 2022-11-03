# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("GenomicRanges")
  library("ArchR")
})

# folders
peaks = "../data/CutNTag/bed/"

# find overlapping peaks between to narrowpeak sets
# extend TRUE means peaks are extended
overlap = function(peak_set1, peak_set2, extend = TRUE) {
  peak_set1 = fread(glue("{peaks}{peak_set1}"))
  peak_set1 = GRanges(
    seqnames = peak_set1$V1,
    ranges = IRanges(start = peak_set1$V2,
                     end = peak_set1$V3)
  )
  
  peak_set2 = fread(glue("{peaks}{peak_set2}"))
  peak_set2 = GRanges(
    seqnames = peak_set2$V1,
    ranges = IRanges(start = peak_set2$V2,
                     end = peak_set2$V3)
  )
  
  if(extend) {
    peak_set2 = extendGR(peak_set2, upstream = 1000, downstream = 1000)
  }
  
  ol = findOverlaps(peak_set1,
                    peak_set2,
                    type = "any",
                    ignore.strand = FALSE)
  ol = only_bulk = peak_set1[queryHits(ol)]
  ol = as.data.frame(ol)
  ol = ol[, 1:3]
  
  return(ol)
  
}

# common G4 and CTCF top peaks (peak score > 75th percentile!)
list.files(peaks)
common_g4_ctcf = overlap(peak_set1 = "G4_AID_0h_quant75_peaks.bed", peak_set2 = "CTCF_AID_0h_quant75_peaks.bed")
write_tsv(common_g4_ctcf, glue("{peaks}common_CTCF_G4_AID_0h_quant75_peaks.bed"), col_names = FALSE)



'
