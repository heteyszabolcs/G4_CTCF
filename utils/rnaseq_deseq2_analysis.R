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
})

# folders
result_folder = "../results/rna_seq_deseq2/"

# prepare tables
counts = fread("../data/RNA-Seq/nextflow_output/star_rsem/rsem.merged.transcript_counts.tsv")
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
res = res %>% drop_na("padj") %>% arrange(padj)

write_tsv(res, glue("{result_folder}RNA-Seq_treat_vs_contr_DESeq2_fc.tsv"))

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
aggr_bed %>% dplyr::filter(group == "high medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(.,
            glue("{result_folder}NT_highmedium.bed"),
            col_names = FALSE)
aggr_bed %>% dplyr::filter(group == "medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(., glue("{result_folder}NT_medium.bed"), col_names = FALSE)
aggr_bed %>% dplyr::filter(group == "low medium") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(.,
            glue("{result_folder}NT_lowmedium.bed"),
            col_names = FALSE)
aggr_bed %>% dplyr::filter(group == "low") %>%
  dplyr::select(seqname, start_position, end_position, gene_symbol, score, strand) %>%
  mutate(strand = ifelse(strand == -1, "-", "+")) %>% 
  write_tsv(., glue("{result_folder}NT_low.bed"), col_names = FALSE)

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
