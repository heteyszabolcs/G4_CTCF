suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("glue")
  library("wigglescout")
  library("GenomicRanges")
  library("ggpubr")
  library("ggrastr")
  library("plyr")
  library("ComplexHeatmap")
  library("EnsDb.Mmusculus.v79")
  library("topGO")
  library("circlize")
})

# folders
peaks = "../data/CutNTag/bed/"
bigwigs = "../data/CutNTag/bw/"
hic_ins = "../data/Hi-C/bed/"
result_folder = "../results/cutntag/"

# TAD boundaries from HiC experiments
hic_beds = list.files(hic_ins, full.names = TRUE, pattern = "noheader.bed")

# RNA-Seq data
deseq2 = "../results/rna_seq_deseq2/reg_log_norm_counts.tsv"
vst = "../results/rna_seq_deseq2/vst_norm_counts.tsv"
deseq2 = fread(deseq2)
vst = fread(vst)
fc = fread("../results/rna_seq_deseq2/RNA-Seq_treat_vs_contr_DESeq2_fc.tsv")

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

ctcf_1 = fread(glue("{peaks}CTCF_AID_0h_downsampled_peaks.narrowPeak"))
#ctcf_2 = fread(glue("{peaks}CTCF_AID_0h_2_norm.narrowPeak"))
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

# add expression to common peaks
G4_CTCF_annot = fread(glue("{peaks}CTCF_G4_WT_common_peaks_annot.tsv"))
G4_CTCF_annot = G4_CTCF_annot %>% left_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  dplyr::rename(TMPyP4_vs_ctrl_fc = log2FoldChange) %>%
  left_join(., vst, by = c("gene_symbol" = "gene_symbol")) %>% 
  dplyr::rename(norm_expr = "NT_1") %>% 
  dplyr::select(-NT_2, -TMPyP4_1, -TMPyP4_2, -lfcSE, -stat, -pvalue, -baseMean) %>% 
  dplyr::filter(abs(DistanceToTSS) < 3000) %>% 
  mutate(peakset = "G4-CTCF") %>% 
  distinct(gene_symbol, .keep_all = TRUE)

CTCF_noG4_annot = fread(glue("{peaks}CTCF_no_G4_WT_peaks_annot.tsv"))
CTCF_noG4_annot = CTCF_noG4_annot %>% left_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  dplyr::rename(TMPyP4_vs_ctrl_fc = log2FoldChange) %>%
  left_join(., vst, by = c("gene_symbol" = "gene_symbol")) %>% 
  dplyr::rename(norm_expr = "NT_1") %>% 
  dplyr::select(-NT_2, -TMPyP4_1, -TMPyP4_2, -lfcSE, -stat, -pvalue, -baseMean) %>% 
  dplyr::filter(abs(DistanceToTSS) < 3000) %>% 
  mutate(peakset = "CTCF only") %>% 
  distinct(gene_symbol, .keep_all = TRUE)

G4_noCTCF = fread(glue("{peaks}G4_no_CTCF_WT_peaks_annot.tsv"))
G4_noCTCF = G4_noCTCF %>% left_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  dplyr::rename(TMPyP4_vs_ctrl_fc = log2FoldChange) %>%
  left_join(., vst, by = c("gene_symbol" = "gene_symbol")) %>% 
  dplyr::rename(norm_expr = "NT_1") %>% 
  dplyr::select(-NT_2, -TMPyP4_1, -TMPyP4_2, -lfcSE, -stat, -pvalue, -baseMean) %>% 
  dplyr::filter(abs(DistanceToTSS) < 3000) %>% 
  mutate(peakset = "G4 only") %>% 
  distinct(gene_symbol, .keep_all = TRUE)

bp = rbind(G4_CTCF_annot, CTCF_noG4_annot, G4_noCTCF) %>% 
  dplyr::select(norm_expr, peakset) %>%
  drop_na()
order = factor(bp$peakset, levels = c("CTCF only", "G4-CTCF", "G4 only"))
meds = ddply(bp, .(peakset), summarise, med = median(norm_expr))

ggplot(bp,
       aes(x = order, y = norm_expr)) +
  geom_boxplot(color = "black", fill = "#fc9272") +
  ylim(0, 20) +
  labs(title = "expression of closest gene (+/- 3 kb to TSS)",
       x = "peak set",
       y = "DESeq2 norm. expr.",
       fill = "") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 25),
    axis.text.x = element_text(
      size = 25,
      color = "black",
      angle = 0
    ),
    axis.text.y = element_text(size = 25, color = "black"),
    axis.title.y = element_text(size = 20,),
    axis.title.x = element_text(size = 20,),
  ) +
  stat_compare_means(label.y = 19,
                     label.x = 1.5,
                     size = 8) +
  stat_compare_means(
    label = "p.signif",
    method = "t.test",
    ref.group = ".all.",
    label.y = 17,
    size = 10
  ) +
  geom_text(data = meds, aes(x = peakset, y = 5, label = round(med, 1)), 
            size = 10, vjust = -1.5)

ggsave(
  plot = last_plot(),
  glue("{result_folder}expr_levels_of_Venn_sets.pdf"),
  width = 8,
  height = 6
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}expr_levels_of_Venn_sets.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# create topGO input. 
go = rbind(G4_CTCF_annot, CTCF_noG4_annot, G4_noCTCF) %>% 
  dplyr::select(TMPyP4_vs_ctrl_fc, gene_symbol, peakset) %>%
  drop_na()
sign_g4_ctcf = go %>% dplyr::filter(peakset == "G4-CTCF") %>% 
  dplyr::filter(abs(TMPyP4_vs_ctrl_fc) > 0.1) %>% pull(gene_symbol)
all_g4_ctcf = go %>% dplyr::filter(peakset == "G4-CTCF") %>% pull(gene_symbol)
sign_g4 = go %>% dplyr::filter(peakset == "G4 only") %>% 
  dplyr::filter(abs(TMPyP4_vs_ctrl_fc) > 0.1) %>% pull(gene_symbol)
all_g4 = go %>% dplyr::filter(peakset == "G4 only") %>% pull(gene_symbol)
sign_ctcf = go %>% dplyr::filter(peakset == "CTCF only") %>% 
  dplyr::filter(abs(TMPyP4_vs_ctrl_fc) > 0.1) %>% pull(gene_symbol)
all_ctcf = go %>% dplyr::filter(peakset == "CTCF only") %>% pull(gene_symbol)

create_input = function(all, sign) {
  input = rep(0, length(all))
  names(input) = all
  input[names(input) %in% sign] = 1
  
  return(input)
}

## topGO analysis
# create topGO output
# with top 10 enriched GO terms
# reference: "org.Mm.eg.db"
create_go_matrix = function(genes, colname) {
  # find biological process ontology
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    # use biological process ontology
    allGenes = genes,
    geneSelectionFun = function(x)
      (x == 1),
    annot = annFUN.org,
    mapping = "org.Mm.eg.db",
    ID = "symbol"
  )
  
  # Fisher test
  resultFisher <-
    runTest(GOdata, algorithm = "elim", statistic = "fisher")
  out <-
    GenTable(GOdata,
             Fisher = resultFisher,
             topNodes = 10,
             numChar = 60)
  
  out = out %>% dplyr::select(Term, Fisher)
  colnames(out) = c("Term", colname)
  
  return(out)
  
}

g4_ctcf_go = create_go_matrix(create_input(all_g4_ctcf, sign_g4_ctcf), colname = "G4-CTCF")
ctcf_only_go = create_go_matrix(create_input(all_ctcf, sign_ctcf), colname = "CTCF only")
g4_only_go = create_go_matrix(create_input(all_g4, sign_g4), colname = "G4 only")

go_outputs = list(ctcf_only_go, g4_ctcf_go, g4_only_go)

full <- Reduce(function(x, y, ...)
  merge(x, y, all = TRUE, ...),
  go_outputs)

full[is.na(full)] = 1
full = full %>% distinct(Term, .keep_all = TRUE) 
rownames(full) = full$Term
full = full %>% dplyr::select(-Term) %>% mutate_if(is.character, as.numeric)
full = as.matrix(full)

col_fun = colorRamp2(c(1, 0.01, 0.001), c("#bdbdbd", "#ffeda0", "#de2d26"))
hm = Heatmap(
  full,
  column_title = "",
  row_title = "",
  name = "Fisher test",
  col = col_fun,
  heatmap_legend_param = list(
    title = "Fisher test",
    at = c(1, -0.05, -0.001),
    labels = c("not enriched", "0.001", "0.01"),
    legend_height = unit(1, "cm")
  ),
  rect_gp = gpar(col = "black", lwd = 0.1),
  show_column_dend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  heatmap_width = unit(7, "cm"),
  heatmap_height = unit(10, "cm"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
hm
pdf(
  file = glue("{result_folder}topGO_on_peaksets.pdf"),
  width = 5,
  height = 5
)
print(hm)
dev.off()

  