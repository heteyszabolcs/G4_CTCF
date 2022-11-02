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

# result table
res = results(dds)
res = as.data.frame(res)
res["gene_name"] = rownames(res)

options(digits = 2)
res = as_tibble(res)
res = res %>% drop_na("padj") %>% arrange(padj)

write_tsv(res, glue("{result_folder}RNA-Seq_treat_vs_contr_DESeq2_fc.tsv"))

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

# volcano plot - thr: p = 0.05, fc = 0.5
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
