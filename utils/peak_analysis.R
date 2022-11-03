suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("glue")
  library("wigglescout")
  library("GenomicRanges")
  library("ggpubr")
  library("ggrastr")
})

# folders
peaks = "../data/CutNTag/bed/"
bigwigs = "../data/CutNTag/bw/"
hic_ins = "../data/Hi-C/bed/"
result_folder = "../results/cutntag/"

# TAD boundaries from HiC experiments
hic_beds = list.files(hic_ins, full.names = TRUE, pattern = "noheader.bed")

# create top peak list
get_tops = function(narrowpeak_set, n = 1000) {
  tops = narrowpeak_set %>% top_n(n = n, wt = V5)
  tops = GRanges(seqnames = tops$V1,
                 ranges = IRanges(start = tops$V2,
                                  end = tops$V3))
  
  return(tops)
  
}

# create top insulating boundaries (negative ins. score, higher insulation)
get_insulating = function(hic_bed, n = 1000) {
  most_ins = hic_bed %>% top_n(n = -n, wt = V4)
  most_ins = GRanges(seqnames = most_ins$V1,
                     ranges = IRanges(start = most_ins$V2,
                                      end = most_ins$V3,
                                      score = most_ins$V4))
  
  return(most_ins)
}

# bigwig scatterplot
create_top_scatter = function(bigwig1 = glue("{bigwigs}CTCF_0h_AID_merge.bw"),
                              bigwig2 = glue("{bigwigs}CTCF_6h_AID_merge.bw"),
                              loci) {
  sc = plot_bw_loci_scatter(bigwig1,
                            bigwig2,
                            loci = loci,
                            verbose = FALSE)
  sc = sc + geom_point(color = "#fc9272") +
    ggtitle(" ") +
    xlim (0, 1500) +
    ylim(0, 1500) +
    theme(
      text = element_text(size = 5),
      plot.title = element_text(size = 3),
      axis.text.x = element_text(size = 5, color = "black"),
      axis.text.y = element_text(size = 5, color = "black")
    )
  
  
  return(print(sc))
}

ctcf_1 = fread(glue("{peaks}CTCF_AID_0h_1_norm.narrowPeak"))
ctcf_2 = fread(glue("{peaks}CTCF_AID_0h_2_norm.narrowPeak"))
top_ctcf = get_tops(ctcf_1)

hic_bed = fread("../data/Hi-C/bed/NT_fastq_merge_bd_200kb_noheader.bed")
top_ins = get_insulating(hic_bed)

narrowpeak_files = list.files(peaks)
bigwig_files = list.files(bigwigs, full.names = TRUE)

# systematic scatterplot analysis - top CTCF peaks
aid_bigwigs = bigwig_files[str_detect(bigwig_files, "AID")]
top_sc_list = list()
counter = 0
for (i in 1:length(aid_bigwigs)) {
  for (j in 1:length(aid_bigwigs)) {
    counter = ifelse(i %in% 1:length(aid_bigwigs), counter + 1, counter)
    plot = create_top_scatter(bigwig1 = aid_bigwigs[i],
                              bigwig2 = aid_bigwigs[j],
                              loci = top_ctcf)
    
    top_sc_list[[as.character(counter)]] = plot
  }
}
ggarrange(plotlist = top_sc_list)

ggsave(
  plot = last_plot(),
  glue("{result_folder}aid_treat_topCTCF_sc.pdf"),
  width = 14,
  height = 12
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}aid_treat_topCTCF_sc.png"),
  width = 14,
  height = 12,
  dpi = 300
)

# systematic scatterplot analysis - top insulating boundaries
top_ins_sc_list = list()
counter = 0
for (i in 1:length(aid_bigwigs)) {
  for (j in 1:length(aid_bigwigs)) {
    counter = ifelse(i %in% 1:length(aid_bigwigs), counter + 1, counter)
    plot = create_top_scatter(bigwig1 = aid_bigwigs[i],
                              bigwig2 = aid_bigwigs[j],
                              loci = top_ins)
    
    top_ins_sc_list[[as.character(counter)]] = plot
  }
}
ggarrange(plotlist = top_ins_sc_list)

ggsave(
  plot = last_plot(),
  glue("{result_folder}aid_treat_top_ins_sc.pdf"),
  width = 14,
  height = 12
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}aid_treat_top_ins_sc.png"),
  width = 14,
  height = 12,
  dpi = 300
)
