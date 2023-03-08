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
  library("ggrastr")
  library("pacman")
  library("ggpubr")
})

# helper function
source("C:/Szabolcs/Karolinska/Data/scripts/promoter_annotation.R")

# folders
result_folder = "../results/rna_seq_deseq2/"

# prepare tables
counts = fread("../data/RNA-Seq/mESC_bulk_RNASeq_TMPyP4_study_CTCF-AID/star_rsem/rsem.merged.transcript_counts.tsv")
counts = counts %>% group_by(gene_id) %>% summarise_all(mean) %>% dplyr::select(-transcript_id)
counts = counts %>%
  mutate_if(is.numeric, round)
counts = as.data.frame(counts)
rownames(counts) = counts$gene_id
counts = counts[-1]

coldata = fread("../data/coldata.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
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
write_tsv(norm_counts, glue("{result_folder}median_of_ratios_counts.tsv"))

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
write_tsv(rlog_counts_t_whole, glue("{result_folder}reg_log_norm_counts-full.tsv"))

# regularized log normalization and export
rlog_counts = rlog(dds, blind = FALSE)
rlog_counts = assay(rlog_counts)
rlog_counts_t = as_tibble(rlog_counts)
rlog_counts_t = rlog_counts_t %>% mutate(gene_symbol = rownames(rlog_counts))
write_tsv(rlog_counts_t, glue("{result_folder}reg_log_norm_counts.tsv"))

# vst normalization and export
vst = vst(dds, blind=FALSE)
vst = assay(vst)
vst = as_tibble(vst)
vst = vst %>% mutate(gene_symbol = rownames(rlog_counts))
write_tsv(vst, glue("{result_folder}vst_norm_counts.tsv"))

# result table
res = results(dds)
res = as.data.frame(res)
res["gene_name"] = rownames(res)

options(digits = 2)
res = as_tibble(res)
res_full = res %>% arrange(padj)
res = res %>% drop_na("padj") %>% arrange(padj)

write_tsv(res, glue("{result_folder}RNA-Seq_treat_vs_contr_DESeq2_fc.tsv"))
write_tsv(res_full, glue("{result_folder}RNA-Seq_treat_vs_contr_DESeq2_fc_full.tsv"))

# thresholds: 0.25 log2FC, adj p-value 0.05
# create bed files from diff. expresed genes
res0.25_p0.05 = res %>% dplyr::filter(abs(log2FoldChange) > 0.25 & padj < 0.05)

options(scipen = 999) # turning off scientific numbers
ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://aug2020.archive.ensembl.org")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
location = getBM(
  attributes = c(
    'mgi_symbol',
    "chromosome_name","start_position", "end_position", "strand"
  ),
  mart = ensembl,
  filters = 'mgi_symbol',
  values = res0.25_p0.05$gene_name
)
bed = location %>% mutate(seqname = paste0("chr", as.character(chromosome_name))) %>% 
  mutate(gene_symbol = ".", score = "0") %>% 
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>% 
  mutate(strand = ifelse(strand == -1, "-", "+"))
write_tsv(bed, glue("{result_folder}RNA-Seq_treat_vs_contr_fc0.25_p0.05.bed"), col_names = FALSE)

# assign promoter regions

# group genes based on rlog norm
aggr = rlog_counts_t %>%
  mutate(NT = rowMeans(dplyr::select(., c("NT_1", "NT_2")))) %>%
  mutate(TMPyP4 = rowMeans(dplyr::select(., c(
    "TMPyP4_1", "TMPyP4_2"
  )))) %>%
  dplyr::select(NT, TMPyP4, gene_symbol) %>% mutate(
    group = case_when(
      NT < quantile(NT, .25) ~ "low",
      NT >= quantile(NT, .25) &
        NT < quantile(NT, .50) ~ "low medium",
      NT >= quantile(NT, .50) &
        NT < quantile(NT, .75) ~ "high medium",
      NT >= quantile(NT, .50) ~ "high"
    )
  )

aggr_with_loc = getBM(
  attributes = c(
    'mgi_symbol',
    "chromosome_name","start_position", "end_position", "strand"
  ),
  mart = ensembl,
  filters = 'mgi_symbol',
  values = aggr$gene_symbol
)
aggr_bed = aggr_with_loc %>% inner_join(., aggr, by = c("mgi_symbol" = "gene_symbol")) %>% 
  mutate(seqname = paste0("chr", as.character(chromosome_name))) %>% 
  mutate(gene_symbol = ".", score = "0") %>% 
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand, group) 
  #mutate(promoter = start_position - 3000) %>% # extend to promoter regions
  #dplyr::select(seqname, start_position, end_position, group)

# write grouped genes into bed files
aggr_bed %>% dplyr::filter(group == "high") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(., glue("{result_folder}NT_high.bed"), col_names = FALSE)

add_upstream = aggr_bed %>% dplyr::filter(group == "high") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  mutate(start_position = start_position - 500) %>% 
  write_tsv(., glue("{result_folder}NT_high_500bp_upstream.bed"), col_names = FALSE)

aggr_bed %>% dplyr::filter(group == "high medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(.,
            glue("{result_folder}NT_highmedium.bed"),
            col_names = FALSE)

add_upstream = aggr_bed %>% dplyr::filter(group == "high medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  mutate(start_position = start_position - 500) %>% 
  write_tsv(., glue("{result_folder}NT_highmedium_500bp_upstream.bed"), col_names = FALSE)

aggr_bed %>% dplyr::filter(group == "low medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(.,
            glue("{result_folder}NT_lowmedium.bed"),
            col_names = FALSE)

add_upstream = aggr_bed %>% dplyr::filter(group == "low medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  mutate(start_position = start_position - 500) %>% 
  write_tsv(., glue("{result_folder}NT_lowmedium_500bp_upstream.bed"), col_names = FALSE)

aggr_bed %>% dplyr::filter(group == "low") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(., glue("{result_folder}NT_low.bed"), col_names = FALSE)

add_upstream = aggr_bed %>% dplyr::filter(group == "low") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  mutate(start_position = start_position - 500) %>% 
  write_tsv(., glue("{result_folder}NT_low_500bp_upstream.bed"), col_names = FALSE)

order = factor(aggr$group, levels = c("high", "high medium", "low medium", "low"))
ggplot(aggr, aes(x = order, y = NT, fill = group)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Set3") +
  ylim(0, 20) +
  labs(
    title = "",
    x = "expression levels",
    y = "DESeq2 norm. expression",
    fill = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  ) 



# volcano plot - thr: p = 0.05, log2FC = 1
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
volc = EnhancedVolcano(
  res,
  lab = res$gene_name,
  labSize = 3.5,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 1",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

# volcano plot - thr: p = 0.05, log2FC = 0.5
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.5.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
volc = EnhancedVolcano(
  res,
  lab = res$gene_name,
  labSize = 3.5,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 0.5",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

# volcano plot - thr: p = 0.05, fc = 1
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
volc = EnhancedVolcano(
  res,
  lab = res$gene_name,
  labSize = 3.5,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 1",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

# volcano plot - thr: p = 0.05, fc = 0.25
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
volc = EnhancedVolcano(
  res,
  lab = res$gene_name,
  labSize = 7,
  axisLabSize = 25,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 0.25,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 0.5",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

pdf(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25.pdf"),
  width = 9,
  height = 5
)
volc = EnhancedVolcano(
  res,
  lab = res$gene_name,
  labSize = 6,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 0.25,
  colAlpha = 0.5,
  drawConnectors = TRUE,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 0.25",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

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
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),       
                     limits = c(-1.5, 1.5)) +
  labs(
    title = "TMPyP4 vs. control",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
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
  geom_text_repel(label = labels, size = 6)
ggplot_volc

ggsave(
  glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25_gg.pdf"),
  width = 5,
  height = 5,
  device = "pdf"
)
ggsave(
  glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25_gg.png"),
  width = 5,
  height = 5,
  dpi = 300
)

# volcano plot by enhancedvolcano package
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25.png"),
  width = 10,
  height = 10,
  units = "in",
  res = 300
)
volc = EnhancedVolcano(
  res,
  lab = res$gene_name,
  labSize = 3.5,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 0.25,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 0.25",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

# same differential gene expr analysis, but the control cells are wild-types (not CTCF-AID condition)
# RNA-Seq made on 2023 Febr 27th
# prepare tables
counts = fread("../data/RNA-Seq/mESC_bulk_RNASeq_TMPyP4_study_20230227_WT/star_rsem/rsem.merged.transcript_counts.tsv")
counts = counts %>% group_by(gene_id) %>% summarise_all(mean) %>% dplyr::select(-transcript_id)
counts = counts %>%
  mutate_if(is.numeric, round)
counts = as.data.frame(counts)
rownames(counts) = counts$gene_id
counts = counts[-1]

coldata = fread("../data/coldata.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
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
norm_counts2 = counts(dds, normalized = TRUE)
genes = rownames(norm_counts2)
norm_counts2 = as_tibble(norm_counts2)
norm_counts2 = norm_counts2 %>% mutate(gene_symbol = genes)
write_tsv(norm_counts2, glue("{result_folder}median_of_ratios_counts-20230227.tsv"))

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
write_tsv(rlog_counts_t_whole, glue("{result_folder}reg_log_norm_counts-full-20230227.tsv"))

# regularized log normalization and export
rlog_counts = rlog(dds, blind = FALSE)
rlog_counts = assay(rlog_counts)
rlog_counts_t = as_tibble(rlog_counts)
rlog_counts_t = rlog_counts_t %>% mutate(gene_symbol = rownames(rlog_counts))
write_tsv(rlog_counts_t, glue("{result_folder}reg_log_norm_counts-20230227.tsv"))

# vst normalization and export
vst = vst(dds, blind=FALSE)
vst = assay(vst)
vst = as_tibble(vst)
vst = vst %>% mutate(gene_symbol = rownames(rlog_counts))
write_tsv(vst, glue("{result_folder}vst_norm_counts-20230227.tsv"))

# result table
res2 = results(dds)
res2 = as.data.frame(res2)
res2["gene_name"] = rownames(res2)

options(digits = 2)
res2 = as_tibble(res2)
res_full2 = res2 %>% arrange(padj)
res2 = res2 %>% drop_na("padj") %>% arrange(padj)

write_tsv(res2, glue("{result_folder}RNA-Seq_treat_vs_contr_DESeq2_fc-20230227.tsv"))
write_tsv(res_full2, glue("{result_folder}RNA-Seq_treat_vs_contr_DESeq2_fc_full-20230227.tsv"))

# thresholds: 0.25 log2FC, adj p-value 0.05
# create bed files from diff. expresed genes
res0.25_p0.05 = res2 %>% dplyr::filter(abs(log2FoldChange) > 0.25 & padj < 0.05)

options(scipen = 999) # turning off scientific numbers
ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://aug2020.archive.ensembl.org")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
location = getBM(
  attributes = c(
    'mgi_symbol',
    "chromosome_name","start_position", "end_position", "strand"
  ),
  mart = ensembl,
  filters = 'mgi_symbol',
  values = res0.25_p0.05$gene_name
)
bed = location %>% mutate(seqname = paste0("chr", as.character(chromosome_name))) %>% 
  mutate(gene_symbol = ".", score = "0") %>% 
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>% 
  mutate(strand = ifelse(strand == -1, "-", "+"))
write_tsv(bed, glue("{result_folder}RNA-Seq_treat_vs_contr_fc0.25_p0.05-20230227.bed"), col_names = FALSE)

# add promoter regions
res0.25_p0.05_proms = assign_promoter_region_mm10(unique(res0.25_p0.05$gene_name))
res0.25_p0.05_proms = res0.25_p0.05_proms %>% distinct(., gene_symbol, .keep_all = TRUE) %>% dplyr::select(seqnames, start, end) %>%  
  write_tsv(., glue("{result_folder}RNA-Seq_treat_vs_contr_fc0.25_p0.05-20230227_proms.bed"), col_names = FALSE)

# group genes based on rlog norm
aggr = rlog_counts_t %>%
  mutate(NT = rowMeans(dplyr::select(., c("NT_1", "NT_2")))) %>%
  mutate(TMPyP4 = rowMeans(dplyr::select(., c(
    "TMPyP4_1", "TMPyP4_2"
  )))) %>%
  dplyr::select(NT, TMPyP4, gene_symbol) %>% mutate(
    group = case_when(
      NT < quantile(NT, .25) ~ "low",
      NT >= quantile(NT, .25) &
        NT < quantile(NT, .50) ~ "low medium",
      NT >= quantile(NT, .50) &
        NT < quantile(NT, .75) ~ "high medium",
      NT >= quantile(NT, .50) ~ "high"
    )
  )

aggr_with_loc = getBM(
  attributes = c(
    'mgi_symbol',
    "chromosome_name","start_position", "end_position", "strand"
  ),
  mart = ensembl,
  filters = 'mgi_symbol',
  values = aggr$gene_symbol
)
aggr_bed = aggr_with_loc %>% inner_join(., aggr, by = c("mgi_symbol" = "gene_symbol")) %>% 
  mutate(seqname = paste0("chr", as.character(chromosome_name))) %>% 
  mutate(gene_symbol = ".", score = "0") %>% 
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand, group) 
#mutate(promoter = start_position - 3000) %>% # extend to promoter regions
#dplyr::select(seqname, start_position, end_position, group)

# write grouped genes into bed files
aggr_bed %>% dplyr::filter(group == "high") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(., glue("{result_folder}NT_high-20230227.bed"), col_names = FALSE)

add_upstream = aggr_bed %>% dplyr::filter(group == "high") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  mutate(start_position = start_position - 500) %>% 
  write_tsv(., glue("{result_folder}NT_high_500bp_upstream-20230227.bed"), col_names = FALSE)

aggr_bed %>% dplyr::filter(group == "high medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(.,
            glue("{result_folder}NT_highmedium-20230227.bed"),
            col_names = FALSE)

add_upstream = aggr_bed %>% dplyr::filter(group == "high medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  mutate(start_position = start_position - 500) %>% 
  write_tsv(., glue("{result_folder}NT_highmedium_500bp_upstream-20230227.bed"), col_names = FALSE)

aggr_bed %>% dplyr::filter(group == "low medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(.,
            glue("{result_folder}NT_lowmedium.bed"),
            col_names = FALSE)

add_upstream = aggr_bed %>% dplyr::filter(group == "low medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  mutate(start_position = start_position - 500) %>% 
  write_tsv(., glue("{result_folder}NT_lowmedium_500bp_upstream-20230227.bed"), col_names = FALSE)

aggr_bed %>% dplyr::filter(group == "low") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(., glue("{result_folder}NT_low.bed"), col_names = FALSE)

add_upstream = aggr_bed %>% dplyr::filter(group == "low") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  mutate(start_position = start_position - 500) %>% 
  write_tsv(., glue("{result_folder}NT_low_500bp_upstream-20230227.bed"), col_names = FALSE)

order = factor(aggr$group, levels = c("high", "high medium", "low medium", "low"))
ggplot(aggr, aes(x = order, y = NT, fill = group)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Set3") +
  ylim(0, 20) +
  labs(
    title = "",
    x = "expression levels",
    y = "DESeq2 norm. expression",
    fill = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  ) 

# volcano plot - thr: p = 0.05, log2FC = 1
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-20230227.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
volc = EnhancedVolcano(
  res2,
  lab = res2$gene_name,
  labSize = 3.5,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 1",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

# volcano plot - thr: p = 0.05, log2FC = 0.5
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.5-20230227.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
volc = EnhancedVolcano(
  res2,
  lab = res2$gene_name,
  labSize = 3.5,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 0.5",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

# volcano plot - thr: p = 0.05, fc = 1
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-20230227.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
volc = EnhancedVolcano(
  res2,
  lab = res2$gene_name,
  labSize = 3.5,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 1",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

# volcano plot - thr: p = 0.05, fc = 0.25
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25-20230227.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
volc = EnhancedVolcano(
  res2,
  lab = res2$gene_name,
  labSize = 7,
  axisLabSize = 25,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 0.25,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 0.5",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

pdf(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25-20230227.pdf"),
  width = 9,
  height = 5
)
volc = EnhancedVolcano(
  res2,
  lab = res2$gene_name,
  labSize = 6,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 0.25,
  colAlpha = 0.5,
  drawConnectors = TRUE,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 0.25",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

# volcano plot by ggplot2
volc_input = res2 %>% mutate(group = case_when(
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
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-1.5, 1.5, 0.5)),       
                     limits = c(-1.5, 1.5)) +
  labs(
    title = "TMPyP4 vs. control (wild-type)",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
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
  geom_text_repel(label = labels, size = 6)
ggplot_volc

ggsave(
  glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25_gg-20230227.pdf"),
  width = 5,
  height = 5,
  device = "pdf"
)
ggsave(
  glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25_gg-20230227.png"),
  width = 5,
  height = 5,
  dpi = 300
)

# volcano plot by enhancedvolcano package
png(
  file = glue("{result_folder}RNA-Seq_treat_vs_contr_volcano-fc0.25-20230227.png"),
  width = 10,
  height = 10,
  units = "in",
  res = 300
)
volc = EnhancedVolcano(
  res2,
  lab = res2$gene_name,
  labSize = 3.5,
  title = "TMPyP4 - control",
  pCutoff = 0.05,
  FCcutoff = 0.25,
  colAlpha = 0.5,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 0.25",
  col = c('grey', 'grey', 'grey', '#f03b20')
)
volc
dev.off()

## comparison of two result table
sign_genes = res %>% dplyr::filter(abs(log2FoldChange) > 0.25 & padj < 0.05) %>% pull(gene_name) %>% unique

res_comb = res %>% inner_join(., res2, by = "gene_name") %>% 
  dplyr::select(gene_name, FC_CTCF_AID = log2FoldChange.x, FC_wildtype = log2FoldChange.y) %>% 
  mutate(sign_label = ifelse(gene_name %in% sign_genes, "diff. in CTCF-AID comp.", "non-diff. in CTCF-AID comp.")) %>% 
  distinct(gene_name, .keep_all = TRUE)

ggplot(res_comb, aes(x = FC_wildtype, y = FC_CTCF_AID)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = sign_label)
  ) +
  ylim(-1.5,1.5) +
  xlim(-1.5, 1.5) +
  scale_fill_manual(values = c("#fc9272", "#f0f0f0")) +
  scale_color_manual(values = "black") +
  labs(
    title = "",
    x = "log2FoldChange (wild-type ctrl.)",
    y = "log2FoldChange (CTCF-AID ctrl.)",
    fill = ""
  ) +
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
  stat_cor(method = "pearson", label.x = -1, label.y = 1.4, size = 6, p.accuracy = 0.001)

ggsave(
  glue("{result_folder}RNA-Seq_FC_wt_Vs._FC_CTCF-AID_1.pdf"),
  width = 7,
  height = 7,
  device = "pdf"
)
ggsave(
  glue("{result_folder}RNA-Seq_FC_wt_Vs._FC_CTCF-AID_1.png"),
  width = 7,
  height = 7,
  dpi = 300
)

ggplot(res_comb %>% dplyr::filter(sign_label == "diff. in CTCF-AID comp."), aes(x = FC_wildtype, y = FC_CTCF_AID)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = sign_label)
  ) +
  ylim(-1.5,1.5) +
  xlim(-1.5, 1.5) +
  scale_fill_manual(values = "#fc9272", labels = "diff. in CTCF-AID comp.") +
  scale_color_manual(values = "black") +
  labs(
    title = "",
    x = "log2FoldChange (wild-type ctrl.)",
    y = "log2FoldChange (CTCF-AID ctrl.)",
    fill = ""
  ) +
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
  stat_cor(method = "pearson", label.x = -1, label.y = 1.4, size = 6, p.accuracy = 0.001)

ggsave(
  glue("{result_folder}RNA-Seq_FC_wt_Vs._FC_CTCF-AID_2.pdf"),
  width = 7,
  height = 7,
  device = "pdf"
)
ggsave(
  glue("{result_folder}RNA-Seq_FC_wt_Vs._FC_CTCF-AID_2.png"),
  width = 7,
  height = 7,
  dpi = 300
)


# significant genes of RNA-Seq where the control was WT mESC
sign_genes2 = res2 %>% dplyr::filter(abs(log2FoldChange) > 0.25 & padj < 0.05) %>% pull(gene_name) %>% unique
sign_genes2 = setdiff(sign_genes2, sign_genes)
print(glue("# of new diff. genes:  {as.character(length(sign_genes2))}"))

ggplot(res_comb %>% dplyr::filter(gene_name %in% sign_genes2), aes(x = FC_wildtype, y = FC_CTCF_AID)) +
  ggrastr::geom_point_rast(
    size = 2,
    shape = 21,
    colour = "black",
    aes(fill = sign_label)
  ) +
  ylim(-1.5,1.5) +
  xlim(-1.5, 1.5) +
  scale_fill_manual(values = c("#f0f0f0", "#fc9272"), labels = c("diff. only in WT comp.", "diff. in both")) +
  scale_color_manual(values = "black") +
  labs(
    title = "",
    x = "log2FoldChange (wild-type ctrl.)",
    y = "log2FoldChange (CTCF-AID ctrl.)",
    fill = ""
  ) +
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
  stat_cor(method = "pearson", label.x = -1, label.y = 1.4, size = 6, p.accuracy = 0.001)

ggsave(
  glue("{result_folder}RNA-Seq_FC_wt_Vs._FC_CTCF-AID_3.pdf"),
  width = 7,
  height = 7,
  device = "pdf"
)
ggsave(
  glue("{result_folder}RNA-Seq_FC_wt_Vs._FC_CTCF-AID_3.png"),
  width = 7,
  height = 7,
  dpi = 300
)

sign_genes2_norm = norm_counts2 %>% dplyr::filter(gene_symbol %in% sign_genes2) %>% 
  inner_join(., norm_counts, by = "gene_symbol", suffix = c(".WT", ".CTCF_AID")) %>% 
  dplyr::select(gene_symbol, everything())
x_order = colnames(sign_genes2_norm)[which(colnames(sign_genes2_norm) != "gene_symbol")]
y_order = sign_genes2_norm %>% mutate(mean = rowMeans(dplyr::select(sign_genes2_norm, -gene_symbol), na.rm = TRUE)) %>% 
  arrange(mean) %>%
  pull(gene_symbol)

sign_genes2_norm = pivot_longer(sign_genes2_norm, cols = "NT_1.WT":"TMPyP4_2.CTCF_AID", names_to = "condition",
                                values_to = "norm_expr")
y_order = factor(sign_genes2_norm$gene_symbol, levels = y_order)
x_order = factor(sign_genes2_norm$condition, levels = x_order)


hm_sign_genes2_norm = ggplot(sign_genes2_norm, aes(x = x_order, y = y_order, fill = log2(norm_expr))) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#a6bddb",
    mid = "#FFFFCC",
    high = "#fc9272",
    midpoint = max(log2(sign_genes2_norm$norm_expr)) / 2,
    limits = c(1, 15)
  ) +
  xlab(label = "") +
  ylab(label = "diff. expr. genes (wild-type ctrl.)") +
  labs(fill = "log2 norm", title = "Normalized expression (DESeq2)") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 8,
      angle = 45,
      hjust = 1,
      vjust = 1.1
    ),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 20)
  ) +
  coord_fixed()
print(hm_sign_genes2_norm)

ggsave(
  glue("{result_folder}RNA-Seq_normexpr_wt_Vs_CTCF-AID.pdf"),
  width = 7,
  height = 12,
  device = "pdf"
)

ggsave(
  glue("{result_folder}RNA-Seq_normexpr_wt_Vs_CTCF-AID.png"),
  width = 7,
  height = 12,
  dpi = 300
)


summary = res %>% mutate(
  status = case_when(
    log2FoldChange > 0.25 & padj < 0.05 ~ "up",
    log2FoldChange < -0.25 & padj < 0.05 ~ "down",
    TRUE ~ ""
  )) %>% dplyr::filter(status != "") %>% 
  group_by(status) %>% count()

bar1 = ggplot(summary, aes(x = status, y = n, fill = status)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = c("#f0f0f0", "#addd8e")) +
  ylab("# of diff. expr. genes") +
  ggtitle("condition: TMPyP4 vs. CTCF-AID") +
  xlab("") +
  ylim(0, 100) +
  guides(fill = "none") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 20,
      angle = 0
    ),
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20)
  )

summary2 = res2 %>% mutate(
  status = case_when(
    log2FoldChange > 0.25 & padj < 0.05 ~ "up",
    log2FoldChange < -0.25 & padj < 0.05 ~ "down",
    TRUE ~ ""
  )) %>% dplyr::filter(status != "") %>% 
  group_by(status) %>% count()

bar2 = ggplot(summary2, aes(x = status, y = n, fill = status)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = c("#f0f0f0", "#addd8e")) +
  ylab("# of diff. expr. genes") +
  ggtitle("condition: TMPyP4 vs. mESC WT") +
  xlab("") +
  ylim(0, 100) +
  guides(fill = "none") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 20,
      angle = 0
    ),
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20)
  )

ggarrange(bar1, bar2)
                       
ggsave(
  glue("{result_folder}RNA-Seq_summary_bars.pdf"),
  width = 12,
  height = 7,
  device = "pdf"
)

ggsave(
  glue("{result_folder}RNA-Seq_summary_bars.png"),
  width = 12,
  height = 7,
  dpi = 300
)




