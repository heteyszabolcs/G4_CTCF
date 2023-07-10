suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("wigglescout")
  library("GenomicRanges")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("ggrepel")
  library("ggpubr")
})

# annotation script
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# result folder
result_folder = "../results/cutntag/"

# folders
bw = "../data/CutNTag_ChIP/bw/"
peak = "../data/CutNTag_ChIP/bed/"


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
               subset = glue("{peak}CTCF_AID_0h_downsampled_peaks.bed"),
               output_file = "CTCF_vs_G4_0H_AID.tsv")

# AID treatment
ctcf_aid_ctcf_increasing = calc_bigwig_fc(bigwig = glue("{bw}CTCF_AID_6h_merge.norm.bw"), 
               bg_bigwig = glue("{bw}CTCF_AID_0h_merge_5million.norm.bw"),
               subset = glue("{peak}CTCF_AID_0h_downsampled_peaks.bed"),
               output_file = "CTCFup_AID_6h_vs_0h_at_CTCF_peaks.tsv")

ctcf_aid_ctcf_decreasing = calc_bigwig_fc(bigwig = glue("{bw}CTCF_AID_0h_merge_5million.norm.bw"), 
                          bg_bigwig = glue("{bw}CTCF_AID_6h_merge.norm.bw"),
                          subset = glue("{peak}CTCF_AID_0h_downsampled_peaks.bed"),
                          output_file = "CTCFdown_AID_6h_vs_0h_at_CTCF_peaks.tsv")


g4_aid_g4_increasing = calc_bigwig_fc(bigwig = glue("{bw}G4_6h_AID_5million.norm.bw"), 
                             bg_bigwig = glue("{bw}G4_0h_AID_5million.norm.bw"),
                             subset = glue("{peak}CTCF_AID_0h_downsampled_peaks.bed"),
                             output_file = "G4up_AID_6h_vs_0h_at_CTCF_peaks.tsv")

g4_aid_g4_decreasing = calc_bigwig_fc(bigwig = glue("{bw}G4_0h_AID_5million.norm.bw"), 
                                      bg_bigwig = glue("{bw}G4_6h_AID_5million.norm.bw"),
                                      subset = glue("{peak}CTCF_AID_0h_downsampled_peaks.bed"),
                                      output_file = "G4down_AID_6h_vs_0h_at_CTCF_peaks.tsv")


# Cut&Tag bigwig analysis
# Cut&Tag result folder
result_folder = "../results/cutntag/"

# folders
bw = "../data/CutNRun/bigwig/"
peak = "../data/CutNTag_ChIP-Seq/bed/" # promoters of interest
rna_seq = "../results/rna_seq_deseq2/RNA-Seq_treat_vs_contr_DESeq2_fc_full-20230227.tsv"
proseq = "../data/PRO-Seq/GSE130691_mESC.featureCounts.txt"
bonev_enh_prom = "../data/Hi-C/Bonev_2017_et_al_tableS3-enh_prom_interactions.tsv"
tmpyp4_hic = "../data/Hi-C/bed/TMPyP4_only_fastq_merge_200kb_annot.tsv"
nt_hic = "../data/Hi-C/bed/NT_only_fastq_merge_200kb_annot.tsv"

# G4 - no CTCF promoters showing higher G4 signal upon 6h TMPyP4
tmpyp4_6h_G4noCTCFproms = calc_bigwig_fc(bigwig = glue("{bw}G4_CnT_E14_TMPyP4_6h_RPGC.bigwig"), 
               bg_bigwig = glue("{bw}G4_CnT_E14_NT_RPGC.bigwig"),
               subset = glue("{peak}G4_no_CTCF_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed"),
               output_file = "TMPyP4_6h_vs_NT-G4_noCTCF_proms.tsv")
tmpyp4_6h_G4noCTCFproms = tmpyp4_6h_G4noCTCFproms %>% mutate(type = "G4 only")

# CTCF - no G4 promoters showing higher G4 signal upon 6h TMPyP4
tmpyp4_6h_CTCFnoG4 = calc_bigwig_fc(bigwig = glue("{bw}G4_CnT_E14_TMPyP4_6h_RPGC.bigwig"), 
                           bg_bigwig = glue("{bw}G4_CnT_E14_NT_RPGC.bigwig"),
                           subset = glue("{peak}CTCF_no_G4_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed"),
                           output_file = "TMPyP4_6h_vs_NT-CTCF_noG4_proms.tsv")
tmpyp4_6h_CTCFnoG4 = tmpyp4_6h_CTCFnoG4 %>% mutate(type = "CTCF only")

# common CTCF - G4 promoters showing higher G4 signal upon 6h TMPyP4
tmpyp4_6h_CTCFG4common = calc_bigwig_fc(bigwig = glue("{bw}G4_CnT_E14_TMPyP4_6h_RPGC.bigwig"), 
                                    bg_bigwig = glue("{bw}G4_CnT_E14_NT_RPGC.bigwig"),
                                    subset = glue("{peak}CTCF_G4_common_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed"),
                                    output_file = "TMPyP4_6h_vs_NT-CTCF_G4_common_proms.tsv")
tmpyp4_6h_CTCFG4common = tmpyp4_6h_CTCFG4common %>% mutate(type = "CTCF & G4")


tmpyp4_6h_CTCFnoG4_proms = tmpyp4_6h_CTCFnoG4 %>% dplyr::filter(abs(distanceToTSS) < 3000)
write_tsv(tmpyp4_6h_CTCFnoG4_proms, glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4_3kb_to_TSS.tsv"))

# visualization of enrichments by bar plot
all = rbind(tmpyp4_6h_G4noCTCFproms, tmpyp4_6h_CTCFnoG4, tmpyp4_6h_CTCFG4common)
bars = all %>% group_by(type, .drop = FALSE) %>% count() %>% 
  ungroup() %>% 
  add_row(type = "CTCF & G4", n = 0) %>% 
  ggplot(., aes(x = reorder(type, -n), y = n)) +
  geom_bar(stat = "identity", color = "black", fill = "#fc9272") +
  labs(
    title = "TMPyP4 G4 Cut&Tag signal quantification",
    x = "promoter set",
    y = "# of gained G4 signal upon TMPyP4",
    fill = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 14, color = "black")
  )
bars

ggsave(
  plot = bars,
  glue("{result_folder}gained_signals_over_proms-bar.pdf"),
  width = 6,
  height = 4
)

ggsave(
  plot = bars,
  glue("{result_folder}gained_signals_over_proms-bar.png"),
  width = 6,
  height = 4,
  dpi = 300
)

# combine with TMPyP4 RNA-Seq data
expr = fread(rna_seq)
proseq_expr = fread(proseq)
proseq_expr = proseq_expr %>% dplyr::select(gene_symbol = "Geneid", starts_with("mESC")) %>% 
  mutate(mean = round(rowMeans(dplyr::select(., starts_with("mESC")), na.rm = TRUE), 2))



tmpyp4_6h_CTCFnoG4_proms = tmpyp4_6h_CTCFnoG4_proms %>% inner_join(., expr, by = c("gene_symbol" = "gene_name"))
sign = tmpyp4_6h_CTCFnoG4_proms %>% dplyr::filter(pvalue < 0.05) %>% dplyr::select(seqnames, start, end, 
                                                                                   CnT_fold_change = fold_change, 
                                                                                   RNA_Seq_fold_change = log2FoldChange,
                                                                                   pvalue, padj, distanceToTSS,
                                                                                   gene_symbol)

ggplot(sign, aes(x = RNA_Seq_fold_change, y = CnT_fold_change, label = gene_symbol)) +
  geom_point(
    size = 6,
    shape = 21,
    colour = "black",
    fill = "#fc9272"
  ) +
  # scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("end", "start")) +
  scale_color_manual(values = "black") +
  labs(
    title = "G4 Cut&Tag signal fc vs. RNA-Seq fc upon 6h TMPyP4",
    x = "RNA-Seq TMPyP4 vs. non-treated",
    y = "TMPyP4 G4 Cut&Tag over CTCF only promoters",
    fill = ""
  ) +
  ylim(0, 10) +
  xlim(-1.5, 1.5) +
  guides(alpha = FALSE, size = FALSE, fill = FALSE, reverse = TRUE) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 15, color = "black")
  ) +
  geom_label_repel() +
  stat_cor(method = "pearson", label.x = -1, label.y = 8.5, size = 5)

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-fc_scatter.pdf"),
  width = 8,
  height = 6
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-fc_scatter.png"),
  width = 8,
  height = 6,
  dpi = 300
)

proseq_expr = proseq_expr %>% inner_join(., tmpyp4_6h_CTCFnoG4, by = "gene_symbol") %>% 
  dplyr::select(gene_symbol, PROSeq_mean = mean, G4_fold_change = fold_change) %>% 
  mutate(label = ifelse(G4_fold_change > 10, gene_symbol, ""))

ggplot(proseq_expr, aes(x = PROSeq_mean, y = G4_fold_change, label = label)) +
  geom_point(
    size = 3,
    shape = 21,
    colour = "black",
    fill = "#fc9272"
  ) +
  # scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("end", "start")) +
  scale_color_manual(values = "black") +
  labs(
    title = "G4 Cut&Tag signal fc vs. PRO-Seq mean",
    x = "PRO-Seq mean read coverage",
    y = "TMPyP4 G4 Cut&Tag over CTCF only promoters",
    fill = ""
  ) +
  ylim(0, 60) +
  xlim(0, 6000) +
  guides(alpha = FALSE, size = FALSE, fill = FALSE, reverse = TRUE) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 15, color = "black")
  ) +
  geom_label_repel() +
  stat_cor(method = "pearson", label.x = 2500, label.y = 50, size = 5)

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-PROSeq_scatter.pdf"),
  width = 8,
  height = 6
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-PROSeq_scatter.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# Bonev et al. enhancer - promoter Hi-C scores
bonev_enh_prom = fread(bonev_enh_prom)

bonev_enh_prom = bonev_enh_prom %>% 
  dplyr::select(gene, hic_score) %>% 
  inner_join(., tmpyp4_6h_CTCFnoG4, by = c("gene" = "gene_symbol")) %>% 
  dplyr::select(gene_symbol = gene, hic_score, G4_fold_change = fold_change) %>% 
  group_by(gene_symbol, G4_fold_change) %>% summarise(mean_hic_score = mean(hic_score)) %>% 
  mutate(label = ifelse(G4_fold_change > 10, gene_symbol, ""))
  
ggplot(bonev_enh_prom, aes(x = mean_hic_score, y = G4_fold_change, label = label)) +
  geom_point(
    size = 3,
    shape = 21,
    colour = "black",
    fill = "#fc9272"
  ) +
  # scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("end", "start")) +
  scale_color_manual(values = "black") +
  labs(
    title = "G4 Cut&Tag signal fc vs. Bonev et al. enh-prom interactions",
    x = "aggregated Hi-C score (Bonev et al.)",
    y = "TMPyP4 G4 Cut&Tag over CTCF only promoters",
    fill = ""
  ) +
  ylim(0, 60) +
  xlim(20, 80) +
  guides(alpha = FALSE, size = FALSE, fill = FALSE, reverse = TRUE) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 15, color = "black")
  ) +
  geom_label_repel() +
  stat_cor(method = "pearson", label.x = 50, label.y = 50, size = 5)

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-HiC_Bonev_scatter.pdf"),
  width = 8,
  height = 6
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-HiC_Bonev_scatter.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# Hi-C data (200 kb bin with insulation scores)
tmpyp4_hic = fread(tmpyp4_hic)
tmpyp4_hic = tmpyp4_hic %>% 
  dplyr::select(-seqnames, -start, -end) %>% 
  inner_join(., tmpyp4_6h_CTCFnoG4, by = c("gene_symbol" = "gene_symbol")) %>% 
  dplyr::select(gene_symbol, log2_insulation_score, G4_fold_change = fold_change, annotation) %>% 
  mutate(label = ifelse(str_detect(annotation, "Promoter"), gene_symbol, "")) %>% 
  mutate(status = ifelse(str_detect(annotation, "Promoter"), "promoter", "non-promoter"))

ggplot(tmpyp4_hic, aes(x = log2_insulation_score, y = G4_fold_change, label = label, fill = status)) +
  geom_point(
    size = 5,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_manual(values = c("#f0f0f0", "#fc9272")) +
  scale_color_manual(values = "black") +
  labs(
    title = "G4 Cut&Tag signal fc vs. TMPyP4 insulation scores",
    x = "TMPyP4 insulation score",
    y = "TMPyP4 G4 Cut&Tag over CTCF only promoters",
    fill = ""
  ) +
  ylim(0, 15) +
  xlim(-1, 1) +
  guides(alpha = FALSE, size = FALSE, reverse = FALSE, label = FALSE) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 15, color = "black")
  ) + geom_label_repel()
  #stat_cor(method = "pearson", label.x = 0.2, label.y = 10, size = 5)

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-HiC_TMPyP4_ins_score.pdf"),
  width = 8,
  height = 6
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-HiC_TMPyP4_ins_score.png"),
  width = 8,
  height = 6,
  dpi = 300
)


nt_hic = fread(nt_hic)
nt_hic = nt_hic %>% 
  dplyr::select(-seqnames, -start, -end) %>% 
  inner_join(., tmpyp4_6h_CTCFnoG4, by = c("gene_symbol" = "gene_symbol")) %>% 
  dplyr::select(gene_symbol, log2_insulation_score, G4_fold_change = fold_change, annotation) %>% 
  mutate(label = ifelse(str_detect(annotation, "Promoter"), gene_symbol, "")) %>% 
  mutate(status = ifelse(str_detect(annotation, "Promoter"), "promoter", "non-promoter"))

ggplot(nt_hic, aes(x = log2_insulation_score, y = G4_fold_change, label = label, fill = status)) +
  geom_point(
    size = 5,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_manual(values = c("#f0f0f0", "#fc9272")) +
  scale_color_manual(values = "black") +
  labs(
    title = "G4 Cut&Tag signal fc vs. non-treated insulation scores",
    x = "insulation score",
    y = "TMPyP4 G4 Cut&Tag over CTCF only promoters",
    fill = ""
  ) +
  ylim(0, 15) +
  xlim(-1, 1) +
  guides(alpha = FALSE, size = FALSE, reverse = FALSE, label = FALSE) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 15, color = "black")
  ) + geom_label_repel()

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-HiC_ins_score.pdf"),
  width = 8,
  height = 6
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}TMPyP4_6h_vs_NT-CTCF_noG4-HiC_ins_score.png"),
  width = 8,
  height = 6,
  dpi = 300
)
