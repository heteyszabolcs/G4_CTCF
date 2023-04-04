# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("ggrepel")
  library("ggpubr")
  library("DGEobj.utils")
  library("ggsignif")
})

# turning off scientific notation
options(scipen=999)

# call annotation function
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# result folder
result_folder = "../results/rna_seq_deseq2/"
# peak folder
cnt_folder = "../data/CutNTag_ChIP-Seq/bed/"
hic_folder = "../data/Hi-C/"
proseq_folder = "../data/PRO-Seq/"
# RNA-Seq data
deseq2 = "../results/rna_seq_deseq2/reg_log_norm_counts-full-20230227.tsv"
fc = fread(glue("{result_folder}RNA-Seq_treat_vs_contr_DESeq2_fc_full-20230227.tsv"))

# PRO-Seq data (GSE130691)
pro = fread(glue("{proseq_folder}GSE130691_mESC.featureCounts.txt"))
pro = pro %>% dplyr::select(-Chr, -Start, -End)
pro_rc = pro %>% dplyr::select(starts_with("mESC")) %>% as.matrix
pro_tpm = as_tibble(convertCounts(pro_rc, unit = "TPM", geneLength = pro$Length, normalize = "none"))
pro_tpm = pro_tpm %>% mutate(gene_symbol = pro$Geneid, strand = pro$Strand) %>% 
  mutate(., TPM_mean = rowMeans(dplyr::select(., starts_with("mESC")), na.rm = TRUE)) %>% 
  dplyr::select(gene_symbol, strand, starts_with("mESC"), TPM_mean) %>% 
  mutate(strand = ifelse(str_detect(pattern = "\\+", string = strand), "+", "-"))

# essential genes (Shahar et al., 2019) - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6836742/
ess = fread("../data/Shahar_et_al_2019-CRISPR_score_mESC_essentiality.tsv")
ess_genes = ess %>% dplyr::filter(`Fold change day 18` <= -2 &
                                    `FDR corrected p-value` < 0.05) %>%
  pull(`Gene name`)

# promoters
CTCF_G4_common = fread(glue("{cnt_folder}CTCF_G4_common_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed"))
G4_no_CTCF = fread(glue("{cnt_folder}G4_no_CTCF_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed"))
CTCF_no_G4 = fread(glue("{cnt_folder}CTCF_no_G4_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed"))

# enhancer-promoter interactions (Bonev 2017)
enh_prom = fread(glue(
  "{hic_folder}Bonev_2017_et_al_tableS3-enh_prom_interactions.tsv"
))
enh_prom = enh_prom %>% dplyr::select(gene, hic_score)

# annotate promoters
annotate_promoters = function(promoter_bed) {
  annot = mm10_annotation(
    regions = promoter_bed,
    seqname_col = "V1",
    start_col = "V2",
    end_col = "V3",
    feature_1 = NULL,
    feature_2 = NULL
  )
  annot = annot %>% dplyr::select(seqnames, start, end, gene = SYMBOL)
  return(annot)
  
}

# annotation of promoter variations
CTCF_G4_common_annot = annotate_promoters(CTCF_G4_common)
CTCF_G4_common_annot = CTCF_G4_common_annot %>% left_join(., enh_prom, by = c("gene" = "gene"))
CTCF_G4_common_enh_prom = CTCF_no_G4_annot %>% inner_join(., enh_prom, by = c("gene" = "gene"))
CTCF_G4_common_enh_prom = unique(CTCF_G4_common_enh_prom$gene)

CTCF_no_G4_annot = annotate_promoters(CTCF_no_G4)
CTCF_no_G4_annot = CTCF_no_G4_annot %>% left_join(., enh_prom, by = c("gene" = "gene"))
CTCF_no_G4_enh_prom = CTCF_no_G4_annot %>% inner_join(., enh_prom, by = c("gene" = "gene"))
CTCF_no_G4_enh_prom = unique(CTCF_no_G4_enh_prom$gene)

G4_no_CTCF_annot = annotate_promoters(G4_no_CTCF)
G4_no_CTCF_annot = G4_no_CTCF_annot %>% left_join(., enh_prom, by = c("gene" = "gene"))
G4_no_CTCF_enh_prom = CTCF_no_G4_annot %>% inner_join(., enh_prom, by = c("gene" = "gene"))
G4_no_CTCF_enh_prom = unique(G4_no_CTCF_enh_prom$gene)

all_gene = unique(c(CTCF_G4_common_annot$gene, CTCF_no_G4_annot$gene, G4_no_CTCF_annot$gene))
print(glue("Number of promoters with G4, CTCF or G4&CTCF: {as.character(length(all_gene))}"))
double_neg_proms = unique(fc$gene_name[which(!fc$gene_name %in% all_gene)])
print(glue("Number of \"empty\" promoters: {as.character(length(double_neg_proms))}"))

double_neg_proms_enh_prom = unique(double_neg_proms[which(double_neg_proms %in% enh_prom$gene)])

# MA plot function
make_ma = function(promoter_bed, title = "", annotation = TRUE) {
  
  if (annotation) {
    annot_prom = annotate_promoters(promoter_bed = promoter_bed)
    labels = unique(annot_prom$gene)
  } else {
    labels = promoter_bed
  }
  
  fc = fc %>% mutate(label = ifelse(gene_name %in% labels, "1", "0"))
  fc = fc %>% mutate(gene_label = ifelse(abs(log2FoldChange) >= 1.25 &
                                           label == "1", gene_name, ""))
  
  ma =
    ggplot(fc %>% arrange(label), aes(log2(baseMean), log2FoldChange)) +
    geom_point(aes(color = label), size = 0.2) +
    #labs(title = title) +
    ylim(-5, 5) +
    scale_color_manual(values = c("#f0f0f0", "#fc9272")) +
    #guides(color = "none") +
    labs(color = title) +
    geom_text_repel(
      label = fc %>% arrange(label) %>% pull(gene_label), max.overlaps = 50, size = 1.75, segment.size = 0.25
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 13, color = "black"),
      axis.text.y = element_text(size = 13, color = "black")
    )
  print(ma)
  
}

# create MA plots
CTCF_G4_common_ma = make_ma(CTCF_G4_common, title = "CTCF & G4 common")
G4_no_CTCF_ma = make_ma(G4_no_CTCF, title = "G4 only")
CTCF_no_G4_ma = make_ma(CTCF_no_G4, title = "CTCF only")
wo_G4_wo_CTCF_ma = make_ma(double_neg_proms, title = "w/o CTCF, w/o G4", annotation = FALSE)

mas = list(CTCF_G4_common_ma, G4_no_CTCF_ma, CTCF_no_G4_ma, wo_G4_wo_CTCF_ma)
mas = ggarrange(plotlist = mas)

# export
ggsave(
  plot = mas,
  glue("{result_folder}MA-CTCF_G4_promoters.pdf"),
  width = 12,
  height = 8,
  device = "pdf"
)

ggsave(
  plot = mas,
  glue("{result_folder}MA-CTCF_G4_promoters.png"),
  width = 12,
  height = 8
)

# MA plot function with Bonev et al. enhancer-promoter highlights
make_ma_enh_prom = function(promoter_bed, title = "", enh_prom_list, annotation = TRUE) {
  
  if (annotation) {
    annot_prom = annotate_promoters(promoter_bed = promoter_bed)
    labels = unique(annot_prom$gene)
  } else {
    labels = promoter_bed
  }
  
  fc = fc %>% mutate(label = ifelse(gene_name %in% labels, "1", ""))
  fc = fc %>% mutate(enh_prom_status = ifelse(gene_name %in% enh_prom_list &
                                                label == "1", "1", "0"))
  fc = fc %>% mutate(gene_label = ifelse(
    abs(log2FoldChange) >= 1.25 & enh_prom_status == "1",
    gene_name,
    ""
  ))
  
  ma =
    ggplot(fc %>% arrange(label), aes(log2(baseMean), log2FoldChange)) +
    geom_point(aes(color = enh_prom_status), size = 0.2) +
    labs(title = title) +
    ylim(-5, 5) +
    scale_color_manual(values = c("#f0f0f0", "#fc9272")) +
    #guides(color = "none") +
    labs(color = paste("enh-promoter interaction \n (Bonev et al.)")) +
    geom_text_repel(
      label = fc %>% arrange(label) %>% pull(gene_label), max.overlaps = 50, size = 1.75, segment.size = 0.25
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 13, color = "black"),
      axis.text.y = element_text(size = 13, color = "black")
    )
  print(ma)
  
}

# promoter variations
CTCF_G4_common_ma_enh = make_ma_enh_prom(promoter_bed = CTCF_G4_common,
                                         title = "CTCF & G4 common",
                                         enh_prom_list = CTCF_G4_common_enh_prom)
G4_no_CTCF_common_ma_enh = make_ma_enh_prom(promoter_bed = G4_no_CTCF,
                                            title = "G4 only",
                                            enh_prom_list = G4_no_CTCF_enh_prom)
CTCF_no_G4_common_ma_enh = make_ma_enh_prom(promoter_bed = CTCF_no_G4,
                                            title = "CTCF only",
                                            enh_prom_list = CTCF_no_G4_enh_prom)
wo_G4_wo_CTCF_ma_enh = make_ma_enh_prom(double_neg_proms, enh_prom_list = double_neg_proms_enh_prom , title = "w/o CTCF, w/o G4", annotation = FALSE)

mas_enh = list(CTCF_G4_common_ma_enh,
               G4_no_CTCF_common_ma_enh,
               CTCF_no_G4_common_ma_enh,
               wo_G4_wo_CTCF_ma_enh)
mas_enh = ggarrange(plotlist = mas_enh)

# export
ggsave(
  plot = mas_enh,
  glue("{result_folder}MA-CTCF_G4_promoters_Bonev_enh_proms.pdf"),
  width = 12,
  height = 12,
  device = "pdf"
)

ggsave(
  plot = mas_enh,
  glue("{result_folder}MA-CTCF_G4_promoters_Bonev_enh_proms.png"),
  width = 12,
  height = 12
)

# boxplot
norm_expr = fread(deseq2)

CTCF_G4_common_annot = annotate_promoters(CTCF_G4_common)
CTCF_G4_common_annot = CTCF_G4_common_annot %>% inner_join(., norm_expr, by = c("gene" = "gene_symbol")) %>%
  distinct(gene, .keep_all = TRUE) %>% mutate(promoter_set = "CTCF & G4") %>% dplyr::select("gene", "NT_1", "NT_2", "TMPyP4_1", "TMPyP4_2", "promoter_set")

CTCF_no_G4_annot = annotate_promoters(CTCF_no_G4)
CTCF_no_G4_annot = CTCF_no_G4_annot %>% inner_join(., norm_expr, by = c("gene" = "gene_symbol")) %>%
  distinct(gene, .keep_all = TRUE) %>% mutate(promoter_set = "CTCF only") %>% dplyr::select("gene", "NT_1", "NT_2", "TMPyP4_1", "TMPyP4_2", "promoter_set")

G4_no_CTCF_annot = annotate_promoters(G4_no_CTCF)
G4_no_CTCF_annot = G4_no_CTCF_annot %>% inner_join(., norm_expr, by = c("gene" = "gene_symbol")) %>%
  distinct(gene, .keep_all = TRUE) %>% mutate(promoter_set = "G4 only") %>% dplyr::select("gene", "NT_1", "NT_2", "TMPyP4_1", "TMPyP4_2", "promoter_set")

double_neg_proms = tibble(genes_without_G4_CTCF = double_neg_proms) %>% inner_join(., norm_expr, by = c("genes_without_G4_CTCF" = "gene_symbol")) %>%
  distinct(genes_without_G4_CTCF, .keep_all = TRUE) %>% mutate(promoter_set = "w/o CTCF, w/o G4") %>% dplyr::rename(gene = genes_without_G4_CTCF)

norm_expr_proms = bind_rows(CTCF_G4_common_annot, CTCF_no_G4_annot, G4_no_CTCF_annot, double_neg_proms)
norm_expr_proms = pivot_longer(
  norm_expr_proms,
  cols = "NT_1":"TMPyP4_2",
  names_to = "condition",
  values_to = "expr"
)

bp = ggplot(norm_expr_proms, aes(x = promoter_set, y = expr, fill = condition)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  ylim(0, 20) +
  labs(title = "DESeq2 normalized expr. of promoter subsets",
       x = "promoter set",
       y = "DESeq2 norm. expr.",
       fill = "sample") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  )
bp

CTCF_G4_common_annot = annotate_promoters(CTCF_G4_common)
CTCF_G4_common_annot = CTCF_G4_common_annot %>% inner_join(., pro_tpm, by = c("gene" = "gene_symbol")) %>%
  distinct(gene, .keep_all = TRUE) %>% mutate(promoter_set = "CTCF & G4") %>% 
  dplyr::select("gene", "strand", "TPM_mean", "promoter_set")

CTCF_no_G4_annot = annotate_promoters(CTCF_no_G4)
CTCF_no_G4_annot = CTCF_no_G4_annot %>% inner_join(., pro_tpm, by = c("gene" = "gene_symbol")) %>%
  distinct(gene, .keep_all = TRUE) %>% mutate(promoter_set = "CTCF only") %>% 
  dplyr::select("gene", "strand", "TPM_mean", "promoter_set")

G4_no_CTCF_annot = annotate_promoters(G4_no_CTCF)
G4_no_CTCF_annot = G4_no_CTCF_annot %>% inner_join(., pro_tpm, by = c("gene" = "gene_symbol")) %>%
  distinct(gene, .keep_all = TRUE) %>% mutate(promoter_set = "G4 only") %>% 
  dplyr::select("gene", "strand", "TPM_mean", "promoter_set")

double_neg_proms = tibble(genes_without_G4_CTCF = double_neg_proms) %>% 
  inner_join(., pro_tpm, by = c("genes_without_G4_CTCF" = "gene_symbol")) %>%
  distinct(genes_without_G4_CTCF, .keep_all = TRUE) %>% mutate(promoter_set = "w/o CTCF, w/o G4") %>% 
  dplyr::select("gene" = "genes_without_G4_CTCF", "strand", "TPM_mean", "promoter_set")

proseq_proms = bind_rows(CTCF_G4_common_annot, CTCF_no_G4_annot, G4_no_CTCF_annot, double_neg_proms)

my_comparisons <- list( c("w/o CTCF, w/o G4", "CTCF only"), c("w/o CTCF, w/o G4", "CTCF & G4"), 
                        c("w/o CTCF, w/o G4", "G4 only") )
pro_bp = ggplot(proseq_proms, aes(x = promoter_set, y = TPM_mean)) +
  geom_boxplot(fill = "#feb24c") +
  ylim(0, 1300) +
  labs(title = "PRO-Seq output of promoter subsets",
       x = "promoter set",
       y = "TPM",
       fill = "sample") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "w/o CTCF, w/o G4", label.y = 1250, size = 9)

pro_bp



# export
ggsave(
  plot = pro_bp,
  glue("{result_folder}PRO-Seq_CTCF_G4_promoters.pdf"),
  width = 8,
  height = 6,
  device = "pdf"
)

ggsave(
  plot = pro_bp,
  glue("{result_folder}PRO-Seq_CTCF_G4_promoters.png"),
  width = 8,
  height = 6
)

stat_df = compare_means(TPM_mean ~ strand, group.by = "promoter_set", data = proseq_proms, p.adjust.method = "holm") %>%
  mutate(y_pos = 20, p.adj = format.pval(p.adj, digits = 2))

strandspec_pro_bp = ggplot(proseq_proms, aes(x = promoter_set, y = TPM_mean, fill = strand)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#fc9272", "#9ecae1")) +
  ylim(0, 1300) +
  labs(title = "PRO-Seq output of promoter subsets",
       x = "promoter set",
       y = "TPM",
       fill = "sample") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  )
strandspec_pro_bp

# export
ggsave(
  plot = strandspec_pro_bp,
  glue("{result_folder}PRO-Seq_strand-CTCF_G4_promoters.pdf"),
  width = 8,
  height = 6,
  device = "pdf"
)

ggsave(
  plot = strandspec_pro_bp,
  glue("{result_folder}PRO-Seq_strand-CTCF_G4_promoters.png"),
  width = 8,
  height = 6
)


tmpyp4_expr_proms = norm_expr_proms %>% dplyr::filter(str_detect(condition, "TMP")) %>% 
  group_by(promoter_set, gene) %>% summarise(expr = mean(expr))

# boxplot of the TMPyP4 expression levels only + statistics
bp2 = ggplot(tmpyp4_expr_proms, aes(x = promoter_set, y = expr, fill = promoter_set)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  ylim(0, 20) +
  labs(title = "TMPyP4 DESeq2 normalized expr. of promoter subsets", subtitle = "(w/o CTCF, w/o G4 as reference)",
       x = "promoter set",
       y = "DESeq2 norm. expr.",
       fill = "sample") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(label.y = 20, label.x = 2.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "w/o CTCF, w/o G4", label.y = 18.5)

bp2

# export
ggsave(
  plot = bp,
  glue("{result_folder}promoter_sets-norm_expressions.pdf"),
  width = 7,
  height = 5,
  device = "pdf"
)

ggsave(
  plot = bp,
  glue("{result_folder}promoter_sets-norm_expressions.png"),
  width = 7,
  height = 5
)

ggsave(
  plot = bp2,
  glue("{result_folder}promoter_sets-TMPyP4_norm_expressions.pdf"),
  width = 7,
  height = 5,
  device = "pdf"
)

ggsave(
  plot = bp2,
  glue("{result_folder}promoter_sets-TMPyP4_norm_expressions.png"),
  width = 7,
  height = 5
)

# barplot: overlap with essential genes (Shahar et al. 2019)
CTCF_G4_common_annot = annotate_promoters(CTCF_G4_common)
CTCF_no_G4_annot = annotate_promoters(CTCF_no_G4)
G4_no_CTCF_annot = annotate_promoters(G4_no_CTCF)

all_gene = unique(c(CTCF_G4_common_annot$gene, CTCF_no_G4_annot$gene, G4_no_CTCF_annot$gene))
double_neg_proms = unique(fc$gene_name[which(!fc$gene_name %in% all_gene)])

ess_tibble = tibble(
  promoter_set = c("CTCF & G4", "CTCF only", "G4 only", "w/o CTCF w/o G4"),
  number_of_essentials = c(length(intersect(
    unique(CTCF_G4_common_annot$gene), ess_genes
  )),
  length(intersect(
    unique(CTCF_no_G4_annot$gene), ess_genes
  )),
  length(intersect(
    unique(G4_no_CTCF_annot$gene), ess_genes
  )),
  length(intersect(
    unique(double_neg_proms), ess_genes
  )))
)

bar = ggplot(data = ess_tibble, aes(
  x = reorder(promoter_set,-number_of_essentials),
  y = number_of_essentials
)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  ylim(0, 1000) +
  labs(title = "Number of mESC essential genes (Shahar et al., 2019)",
       x = "promoter set",
       y = "") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 16, color = "black", angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 16, color = "black")
  )
bar

# export
ggsave(
  plot = bar,
  glue(
    "{result_folder}promoter_sets-Shahar_etal_mESC_essential_genes.pdf"
  ),
  width = 7,
  height = 5,
  device = "pdf"
)

ggsave(
  plot = bar,
  glue(
    "{result_folder}promoter_sets-Shahar_etal_mESC_essential_genes.png"
  ),
  width = 7,
  height = 5
)
