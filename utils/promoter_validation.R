# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("GenomicRanges")
})

source("C:/Szabolcs/Karolinska/Data/scripts/promoter_annotation.R")

txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
genes = transcriptsBy(txdb, "gene")
proms = promoters(genes, upstream = 1500, downstream = 500)
proms = as_tibble(proms)
proms = proms %>% distinct(seqnames, start, end, .keep_all = TRUE)

mm10_proms = GRanges(
  seqnames = proms$seqnames,
  ranges = IRanges(
    start = proms$start,
    end = proms$end,
    names = rep("promoter - TxDb.Mmusculus.UCSC.mm10.knownGene", length(proms$seqnames)),
  )
)

# our merged promoter set
merged_proms = "../data/promoters_geneSymbol.intersect.bed"
merged_proms = fread(merged_proms)

g4_only_proms = "../data/CutNTag_ChIP-Seq/bed/G4_no_CTCF_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed"
g4_only_proms = fread(g4_only_proms)
ctcf_only_proms = "../data/CutNTag_ChIP-Seq/bed/CTCF_no_G4_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed"
ctcf_only_proms = fread(ctcf_only_proms)

# merged overlap
merged_proms = GRanges(
  seqnames = merged_proms$V1,
  ranges = IRanges(
    start = merged_proms$V2,
    end = merged_proms$V3,
    names = rep("merged", length(merged_proms$V1)),
  )
)

merged_proms_ol = findOverlaps(merged_proms, mm10_proms, minoverlap = 1,
                          type = "any",
                          ignore.strand = FALSE)

overlapping_proms = as_tibble(merged_proms[queryHits(merged_proms_ol)])
overlapping_proms = overlapping_proms %>% distinct(seqnames, start, end, .keep_all = TRUE) %>% 
  dplyr::select(seqnames, start, end) %>% 
  write_tsv("../results/cutntag/promoters_overlap_with_mm10_set.tsv", col_names = FALSE)

merged_proms_ol = as.data.frame(merged_proms_ol)
merged_proms_ol = unique(merged_proms_ol$queryHits)

merged_proms_perc = round((length(merged_proms_ol) / length(merged_proms)) * 100, 1)
merged_proms = as_tibble(merged_proms)
length(rownames(merged_proms)) - length(merged_proms_ol)
length(merged_proms_ol) - length(mm10_proms)

# G4 only overlap
g4_only_proms = GRanges(
  seqnames = g4_only_proms$V1,
  ranges = IRanges(
    start = g4_only_proms$V2,
    end = g4_only_proms$V3,
    names = rep("G4 only promoters", length(g4_only_proms$V1)),
  )
)

g4_only_ol = findOverlaps(g4_only_proms, mm10_proms, minoverlap = 1,
             type = "any",
             ignore.strand = FALSE)
g4_only_ol = as.data.frame(g4_only_ol)
g4_only_ol = unique(g4_only_ol$queryHits)

g4_only_perc = round((length(g4_only_ol) / length(g4_only_proms)) * 100, 1)
g4_only_proms = as_tibble(g4_only_proms)
length(rownames(g4_only_proms)) - length(g4_only_ol)
length(mm10_proms) - length(g4_only_ol)

# CTCF only overlap
ctcf_only_proms = GRanges(
  seqnames = ctcf_only_proms$V1,
  ranges = IRanges(
    start = ctcf_only_proms$V2,
    end = ctcf_only_proms$V3,
    names = rep("CTCF only promoters", length(ctcf_only_proms$V1)),
  )
)

ctcf_only_ol = findOverlaps(ctcf_only_proms, mm10_proms, minoverlap = 1,
                          type = "any",
                          ignore.strand = FALSE)
ctcf_only_ol = as.data.frame(ctcf_only_ol)
ctcf_only_ol = unique(ctcf_only_ol$queryHits)

ctcf_only_perc = round((length(ctcf_only_ol) / length(ctcf_only_proms)) * 100, 1)
ctcf_only_proms = as_tibble(ctcf_only_proms)
length(rownames(ctcf_only_proms)) - length(ctcf_only_ol)
length(mm10_proms) - length(ctcf_only_ol)



