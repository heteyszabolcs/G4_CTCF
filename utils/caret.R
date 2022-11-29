# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("ggrepel")
  library("ggpubr")
  library("glue")
  library("caret")
  library("caTools")
  library("visreg")
  library("ComplexHeatmap")
  library("circlize")
  
})

set.seed(42)

# result folder
result_folder = "../results/hi-c/"
cnt_folder = "../results/cutntag/"
hic_folder = "../results/hi-c/"
rnaseq_folder = "../results/rna_seq_deseq2/"

# RNA-Seq data
deseq2 = "../results/rna_seq_deseq2/reg_log_norm_counts.tsv"
vst = "../results/rna_seq_deseq2/vst_norm_counts.tsv"
deseq2 = fread(deseq2)
vst = fread(vst)
fc = fread("../results/rna_seq_deseq2/RNA-Seq_treat_vs_contr_DESeq2_fc.tsv")

# Cut&Tag data
g4up = fread(glue("{cnt_folder}G4up_AID_6h_vs_0h_at_CTCF_peaks.tsv"))
g4up = g4up %>% dplyr::rename(G4_foldchange = fold_change) %>% dplyr::rename(G4_distanceToTSS = distanceToTSS)
ctcfdown = fread(glue("{cnt_folder}CTCFdown_AID_6h_vs_0h_at_CTCF_peaks.tsv"))
ctcfdown = ctcfdown %>% dplyr::rename(CTCF_foldchange = fold_change) %>% dplyr::rename(CTCF_distanceToTSS = distanceToTSS)

# find CTCF peaks where G4 enrichment happened with CTCF dropping
# then visualize expression features of genes close to these peaks
g4up_ctcfdown = g4up %>% inner_join(., ctcfdown, by = "gene_symbol") %>% 
  dplyr::filter(abs(G4_foldchange) > 2 & abs(CTCF_foldchange) > 2) %>% 
  dplyr::filter(abs(G4_distanceToTSS) < 3000 &
                  abs(CTCF_distanceToTSS) < 3000) %>% 
  dplyr::mutate(difference = G4_foldchange + CTCF_foldchange) %>% 
  inner_join(., vst) %>% dplyr::select(-NT_2, -TMPyP4_1, -TMPyP4_2) %>% 
  dplyr::rename(norm_expr = NT_1) %>% arrange(desc(norm_expr)) %>% 
  inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>% 
  dplyr::rename(TMPyP4_vs_ctrl_fc = log2FoldChange) %>% 
  arrange(desc(TMPyP4_vs_ctrl_fc))

expr_hm = as.matrix(g4up_ctcfdown$norm_expr)
col_fun = colorRamp2(c(10, 13, 16), c("#ffeda0", "#feb24c", "#f03b20"))
expr_hm = Heatmap(
  expr_hm,
  column_title = "",
  row_title = "closest gene",
  name = "norm. expr.",
  col = col_fun,
  show_column_dend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  heatmap_width = unit(2, "cm"),
  heatmap_height = unit(14, "cm"),
  row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
fc_hm = as.matrix(g4up_ctcfdown$TMPyP4_vs_ctrl_fc)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("#636363", "#f0f0f0", "#feb24c"))
fc_hm = Heatmap(
  fc_hm,
  column_title = "",
  row_title = "closest gene",
  name = "RNA-Seq log2FC",
  col = col_fun,
  show_column_dend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  heatmap_width = unit(2, "cm"),
  heatmap_height = unit(14, "cm"),
  row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
pdf(
  file = glue("{rnaseq_folder}G4up_2fold_CTCFdown_2fold_proms_expr_hm.pdf"),
  width = 5,
  height = 6
)
hm_list = fc_hm + expr_hm
print(hm_list)
dev.off()

# write to bed in order to vis by deeptools
bed = g4up_ctcfdown %>% dplyr::select(seqnames.x, start.x, end.x) %>% write_tsv(., glue("{cnt_folder}G4up_2fold_CTCFdown_2fold_prom_regions.bed"), col_names = FALSE)

# HiC data
tmpyp4_only_board = fread(glue("{hic_folder}TMPyP4_only_fastq_merge_200kb_annot.tsv"))

# random forest
input = g4up %>% inner_join(., vst, by = "gene_symbol") %>%
  inner_join(., tmpyp4_only_board, by = "gene_symbol") %>%
  inner_join(., fc, by = c("gene_symbol" = "gene_name")) %>%
  inner_join(., ctcfdown, by = "gene_symbol") %>%
  dplyr::select(
    gene_symbol,
    G4_distanceToTSS,
    CTCF_distanceToTSS,
    G4_foldchange = fold_change,
    norm_expr_NT = NT_1,
    TMPyP4_log2_ins_score = log2_insulation_score,
    RNA_Seq_log2_foldchange = log2FoldChange,
    CTCF_foldchange
  ) %>%
  mutate(target = case_when(
    G4_foldchange > 2 &
      CTCF_foldchange > 1 & abs(G4_distanceToTSS) < 3000 &
      abs(CTCF_distanceToTSS) < 3000 ~ "1",
    TRUE ~ "0"
  )) %>%
  dplyr::select(-ends_with("distanceToTSS"),
                -G4_foldchange,
                -CTCF_foldchange,
                -gene_symbol)

input = input %>% distinct_all()

control = trainControl(
  method = 'repeatedcv',
  number = 10,
  repeats = 3,
  search = 'random'
)

rf_random = train(
  target ~ .,
  data = input,
  method = 'rf',
  metric = 'Accuracy',
  tuneLength  = 15,
  trControl = control
)

rf_imp = varImp(rf_random, scale = TRUE)

# logistic regression
lr = train(target ~ ., data = input, method = "glm", family = "binomial")
summary(lr)


