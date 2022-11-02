# packages
suppressPackageStartupMessages({
  library("dplyr")
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("ggpubr")
  library("cowplot")
  library("glue")
  library("GenomicRanges")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("Mus.musculus")
})

# add result folder
result_folder = "../results/hi-c/"

# bed files related to TAD boundaries
bed = "../data/Hi-C/bed/"
beds = list.files(bed, full.names = TRUE)
# bed file names (aka experiments)
names = c("common", "NT", "NT_only", "TMPyP4", "TMPyP4_only")

# generate object about mm10 known gene symbols
TxDb(Mus.musculus) = TxDb.Mmusculus.UCSC.mm10.knownGene
mm10_genes = transcriptsBy(Mus.musculus, columns = "SYMBOL")


bed_tibbles = list()
for (i in 1:length(beds)) {
  x = fread(beds[i])
  x = x %>% mutate(treatment = names[i])
  bed_tibbles[[i]] = x
}
bed_tibbles = bind_rows(bed_tibbles)

# check distributions of insulation scores
hist = ggplot(bed_tibbles,
              aes(x = log2_insulation_score_200000, fill = treatment)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(title = "Insulation score distributions",
       x = "log2 insulation score",
       y = "",
       fill = "treatment") +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 10, color = "black")
  )
hist

medians = bed_tibbles %>% group_by(treatment) %>% summarise(median = median(log2_insulation_score_200000)) %>% arrange(desc(median))
order = factor(bed_tibbles$treatment, levels = medians$treatment)

bp = ggplot(bed_tibbles,
            aes(x = order, y = log2_insulation_score_200000, fill = treatment)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  ylim(-10, 4) +
  labs(title = "",
       x = "treatment",
       y = "log2 insulation score",
       fill = "treatment") +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      color = "black",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(size = 10, color = "black")
  ) +
  stat_compare_means(label.y = 3, label.x = 2.5) +
  stat_compare_means(
    label = "p.signif",
    method = "t.test",
    ref.group = ".all.",
    label.y = 2
  )
bp

grid = plot_grid(hist, bp)

ggsave(
  glue("{result_folder}Insulation_score_distr.png"),
  plot = grid,
  width = 8,
  height = 4,
  dpi = 300,
)

# get_intratad function: retrieve TAD domain ranges based on TAD boundary bed files
get_intratad = function(bed_file) {
  bed = fread(bed_file)
  chrs = unique(bed$seqnames)
  
  # retrieve intratad regions
  intratad_ranges = list()
  for (chr in chrs) {
    x_filt = bed %>% dplyr::filter(seqnames == chr)
    tad_start = numeric()
    tad_end = numeric()
    for (i in 1:nrow(x_filt)) {
      tad_start = c(tad_start, x_filt[i, 3]) # TAD start: end of actual region
      tad_end = c(tad_end, x_filt[i + 1, 2]) # TAD end: start of next region
      tibble = tibble(
        seqnames = chr,
        tad_start = as.numeric(tad_start),
        tad_end = as.numeric(tad_end)
      )
    }
    intratad_ranges[[chr]] = tibble
  }
  
  intratad_ranges = bind_rows(intratad_ranges)
  
  # close out NAs (where TAD end is on different chromosome)
  intratad_ranges = intratad_ranges %>%
    na.omit() %>% dplyr::select(seqnames, end = tad_end, start = tad_start)
  return(intratad_ranges)
}

# annotate_intratad function: assign mm10 gene symbols to TAD domains
annotate_intratad = function(tad_regions) {
  tad_regions = as.data.frame(tad_regions)
  tad_regions = makeGRangesFromDataFrame(tad_regions)
  
  # overlapping with known mm10 genes
  hits = findOverlaps(query = mm10_genes,
                      subject = tad_regions,
                      type = "within")
  
  hits = as.data.frame(hits)
  
  seqnames = character()
  starts = character()
  ends = character()
  gene_symbols = character()
  for (i in 1:length(hits$queryHits)) {
    index = hits[i, 1]
    tad_index = hits[i, 2]
    
    gene_symbol = unlist(unique(as.data.frame(mm10_genes[index])$SYMBOL))
    seqname = as.data.frame(tad_regions[tad_index])$seqnames
    start = as.data.frame(tad_regions[tad_index])$start
    end = as.data.frame(tad_regions[tad_index])$end
    
    seqnames = c(seqnames, seqname)
    starts = c(starts, start)
    ends = c(ends, end)
    gene_symbols = c(gene_symbols, gene_symbol)
  }
  
  annotation = tibble(seqnames = seqnames, starts = starts, ends = ends, gene_symbol = gene_symbols)

  return(annotation)
}

# loop through various bed files with boundary informations
common_tad = get_intratad(beds[1])
nt_tad = get_intratad(beds[2])
nt_only_tad = get_intratad(beds[3])
tmpyp4_tad = get_intratad(beds[4])
tmpyp4_only_tad = get_intratad(beds[5])

# annotate them
print("Working on common_200kb_TAD")
common_tad = annotate_intratad(common_tad)
write_tsv(common_tad, glue("{result_folder}common_200kb_TAD_annots.tsv"))

print("Working on NT_200kb_TAD")
nt_tad = annotate_intratad(nt_tad)
write_tsv(nt_tad, "{result_folder}NT_200kb_TAD_annots.tsv")

print("Working on NT_only_200kb_TAD")
nt_only_tad = annotate_intratad(nt_only_tad)
write_tsv(nt_only_tad, glue("{result_folder}NT_only_200kb_TAD_annots.tsv"))

print("Working on TMPyP4_200kb_TAD")
tmpyp4_tad = annotate_intratad(tmpyp4_tad)
write_tsv(tmpyp4_tad, glue("{result_folder}TMPyP4_200kb_TAD_annots.tsv"))

print("Working on TMPyP4_only_200kb")
tmpyp4_only_tad = annotate_intratad(tmpyp4_only_tad)
write_tsv(tmpyp4_only_tad, glue("{result_folder}TMPyP4_only_200kb_TAD_annots.tsv"))


