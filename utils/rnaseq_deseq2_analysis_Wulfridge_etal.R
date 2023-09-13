# packages
suppressPackageStartupMessages({
  library("DESeq2")
  library("data.table")
  library("tidyverse")
  library("EnhancedVolcano")
  library("glue")
  library("enrichR")
  library("scales")
  library("biomaRt")
  library("ggrepel")
  library("ggbreak")
  library("ggrastr")
  library("pacman")
  library("ggpubr")
  library("biomaRt")
})

# helper function
source("C:/Szabolcs/Karolinska/Data/scripts/promoter_annotation.R")

# folders
result_folder = "../results/rna_seq_deseq2/"

# prepare tables
counts = fread("../data/RNA-Seq/Wulfridge_GSE208145/star_rsem/rsem.merged.transcript_counts.tsv")
counts = counts %>% group_by(gene_id) %>% summarise_all(mean) %>% dplyr::select(-transcript_id)
counts = counts %>%
  mutate_if(is.numeric, round)
counts = as.data.frame(counts)
rownames(counts) = counts$gene_id
counts = counts[-1]

coldata = fread("../data/coldata_Wulfridge_RNA-Seq.txt", 
                sep = "\t", 
                header = TRUE, stringsAsFactors = TRUE)
coldata = as.data.frame(coldata)
rownames(coldata) = coldata$sample
all(rownames(coldata) == colnames(counts))

## DESeq2 protocol + differential analysis
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~ treatment)
keep = rowSums(counts(dds)) >= 10
dds = dds[keep, ]
dds$condition = factor(dds$treatment, levels = c("control", "treated"))
dds$condition = relevel(dds$treatment, ref = "control")
dds = DESeq(dds)

# normalized read counts (median of ratios)
norm_counts = counts(dds, normalized = TRUE)
genes = rownames(norm_counts)
norm_counts = as_tibble(norm_counts)
norm_counts = norm_counts %>% mutate(gene_symbol = genes)

# normalized read counts (with low expressed genes)
dds_whole = DESeqDataSetFromMatrix(countData = counts,
                                   colData = coldata,
                                   design = ~ treatment)
dds_whole$condition = factor(dds_whole$treatment, levels = c("control", "treated"))
dds_whole$condition = relevel(dds_whole$treatment, ref = "control")
dds_whole = DESeq(dds_whole)

rlog_counts_whole = rlog(dds_whole, blind = FALSE)
rlog_counts_whole = assay(rlog_counts_whole)
rlog_counts_t_whole = as_tibble(rlog_counts_whole)
rlog_counts_t_whole = rlog_counts_t_whole %>% mutate(gene_symbol = rownames(rlog_counts_whole))

# regularized log normalization and export
rlog_counts = rlog(dds, blind = FALSE)
rlog_counts = assay(rlog_counts)
rlog_counts_t = as_tibble(rlog_counts)
rlog_counts_t = rlog_counts_t %>% mutate(gene_symbol = rownames(rlog_counts))

# vst normalization and export
vst = vst(dds, blind=FALSE)
vst = assay(vst)
vst = as_tibble(vst)
vst = vst %>% mutate(gene_symbol = rownames(rlog_counts))

# result table
res = results(dds)
res = as.data.frame(res)
res["gene_name"] = rownames(res)

options(digits = 2)
res = as_tibble(res)
res_full = res %>% arrange(padj)
res = res %>% drop_na("padj") %>% arrange(padj)

# volcano plot by ggplot2
volc_input = res %>% mutate(group = case_when(
  log2FoldChange > 0.25 & padj < 0.05 ~ "up",
  log2FoldChange < -0.25 & padj < 0.05 ~ "down",
  log2FoldChange >= -0.25 & log2FoldChange <= 0.25 ~ "unaltered"
)) %>% 
  mutate(sign_label = case_when(log2FoldChange > 0.25 & padj < 1e-7 ~ gene_name,
                                log2FoldChange < -0.25 & padj < 1e-7 ~ gene_name,
                                log2FoldChange >= -0.25 & log2FoldChange <= 0.25 ~ ""))
labels = volc_input %>% pull(sign_label)

cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey") 
sizes = c("up" = 4, "down" = 4, "unaltered" = 2) 
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

ggplot_volc = volc_input %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = group,    
             size = group,
             alpha = group)) + 
  ggrastr::geom_point_rast(shape = 21) + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-5, 5, 1)),       
                     limits = c(-5, 5)) +
  ylim(0, 150) +
  scale_y_break(c(40, 120), ticklabels = c(120, 130, 140)) +
  labs(
    title = "Delta Praf2 G4 vs. control (Wulfridge et al.)",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  guides(alpha = "none", size = "none", fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text(label = labels, size = 6, nudge_y = 5)
ggplot_volc

ggsave(
  glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-Wulfridge_et_al.pdf"),
  width = 8,
  height = 5,
  device = "pdf"
)

fc_cutoff = 1
padj_cutoff = 0.05
ups = res %>% dplyr::filter(padj < padj_cutoff & log2FoldChange > fc_cutoff) %>% 
  pull(gene_name) %>% length
downs = res %>% dplyr::filter(padj < padj_cutoff & log2FoldChange < -fc_cutoff) %>% 
  pull(gene_name) %>% length

glue("Number of up-regulated genes: {as.character(ups)}")
glue("Number of down-regulated genes: {as.character(downs)}")


volc_input_praf2 = volc_input %>% mutate(sign_label = ifelse(gene_name == "Praf2", "Praf2", ""))
ggplot_volc_praf2 = volc_input_praf2 %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) + 
  ggrastr::geom_point_rast(shape = 20) + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(-5, 5, 1)),       
                     limits = c(-5, 5)) +
  ylim(0, 25) +
  labs(
    title = "Delta Praf2 G4 vs. control (Wulfridge et al.)",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  guides(alpha = "none", size = "none", fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = volc_input_praf2$sign_label, size = 6, nudge_y = 1, color = "red")
ggplot_volc_praf2

ggsave(
  glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-Wulfridge_et_al_Praf2.pdf"),
  width = 8,
  height = 5,
  device = "pdf"
)

diff_genes = res %>% dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  dplyr::select(gene_name) %>% distinct_all()
ensembl =                                ## fully specified mart
  useMart("ensembl", dataset = "mmusculus_gene_ensembl")  
annot = getBM(
  attributes =  c("mgi_symbol", "chromosome_name"),
  filters =  "mgi_symbol",
  values =  diff_genes$gene_name,
  mart = ensembl
)
res_annot = res %>% left_join(., annot, by = c("gene_name" = "mgi_symbol")) %>% 
  dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  relocate(c("gene_name", "chromosome_name"), .before = baseMean) %>% 
  arrange(padj)

write_tsv(res_annot, 
          glue("{result_folder}Wulfridge_RNA-Seq_Praf2DeltaG4_vs_contr_DESeq2_fc-log2FC1_padj0.05.tsv"))
write_tsv(res, 
          glue("{result_folder}Wulfridge_RNA-Seq_Praf2DeltaG4_vs_contr_DESeq2_fc-full.tsv"))

library("topGO")
create_input = function(final_res, res) {

  genes = final_res %>% 
    pull(gene_name) 
  
  input = rep(0, dim(res)[1])
  names(input) = res$gene_name
  input[names(input) %in% genes] = 1
  
  return(input)
}

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
             topNodes = 20,
             numChar = 60)
  
  out = out %>% dplyr::select(Term, Fisher)
  colnames(out) = c("Term", colname)
  
  return(out)
  
}

input = create_input(final_res = res_annot, res = res)
topgo = create_go_matrix(genes = input, colname = "Fisher_p-value")
write_tsv(topgo, 
          glue("{result_folder}Wulfridge_RNA-Seq_Praf2DeltaG4_vs_contr_DESeq2_fc-topGO.tsv"))
