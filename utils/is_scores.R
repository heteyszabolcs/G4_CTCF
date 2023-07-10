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
  library("wigglescout")
})

# result folder
result_folder = "../results/hi-c/"

# turning off scientific notation
options(scipen = 999)

# call annotation function
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# all insulation score bigwigs (coming from 20 kb Hi-C experiments)
is_bws = c(
  "../data/Hi-C/Nora_AID_NT_all_insulation_score_20kb_win200kb.bw",
  "../data/Hi-C/Nora_AID_IAA_all_insulation_score_20kb_win200kb.bw",
  "../data/Hi-C/Nora_AID_washoff_all_insulation_score_20kb_win200kb.bw",
  "../data/Hi-C/Bonev_insulation_score_20kb_win200kb.bw",
  "../data/Hi-C/NT_fastq_merge_insulation.bw",
  "../data/Hi-C/TMPyP4_fastq_merge_insulation.bw"
)

# 2 kb level bins
is_bins = bw_bins(is_bws,
                  bin_size = 2000,
                  remove_top = 0.1,
                  genome = "mm10",)

is_bins = as_tibble(is_bins)

violins = is_bins %>% dplyr::select(
  "Nora et al. - NT" = "Nora_AID_NT_all_insulation_score_20kb_win200kb",
  "Nora et al. - IAA" = "Nora_AID_IAA_all_insulation_score_20kb_win200kb",
  "Nora et al. - washoff" = "Nora_AID_washoff_all_insulation_score_20kb_win200kb",
  "Bonev et al." = "Bonev_insulation_score_20kb_win200kb",
  "NT" = "NT_fastq_merge_insulation",
  "TMPyP4" = "TMPyP4_fastq_merge_insulation"
) %>% pivot_longer(.,
                   cols = "Nora et al. - NT":"TMPyP4",
                   names_to = "HiC_data",
                   values_to = "aggr_IS")


my_comparisons <-
  list(
    c("TMPyP4", "Nora et al. - IAA"),
    c("TMPyP4", "Nora et al. - NT"),
    c("TMPyP4", "Nora et al. - washoff"),
    c("TMPyP4", "NT"),
    c("TMPyP4", "Bonev et al.")
  )

ggplot(violins, aes(x = HiC_data, y = aggr_IS)) +
  geom_boxplot(color = "black", fill = "#3182bd") +
  ylim(-10, 10) +
  labs(title = "Insulation score distribution, 2 kb bins",
       x = "",
       y = "Aggr. insulation score") +
  guides(fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(
      size = 10,
      color = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 10, color = "black")
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    label.y = c(5, 6, 7, 8, 9, 10),
    method = "t.test",
    label = "p.signif"
  )

# save
ggsave(
  glue(
    "{result_folder}Insulation_scores_200kb_windows_20kb_HiC_res.pdf"
  ),
  plot = last_plot(),
  width = 8,
  height = 6
)

ggsave(
  glue(
    "{result_folder}Insulation_scores_200kb_windows_20kb_HiC_res.png"
  ),
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 500
)

violins_filt = violins %>% 
  dplyr::filter(aggr_IS < quantile(violins$aggr_IS, 0.75)) %>% 
  dplyr::filter(aggr_IS > quantile(violins$aggr_IS, 0.25)) 

order = violins %>% group_by(HiC_data) %>% summarise(median = median(aggr_IS)) %>% 
  arrange(desc(median)) %>% pull(HiC_data) 
order = factor(violins_filt$HiC_data, levels = order)

ggplot(violins_filt, aes(x = order, y = aggr_IS)) +
  geom_boxplot(color = "black", fill = "#3182bd") +
  ylim(-0.30, 0.30) +
  labs(title = "Insulation score distribution, 2 kb bins",
       x = "",
       y = "aggr. insulation score") +
  guides(fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(
      size = 10,
      color = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 10, color = "black")
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    label.y = c(0.12, 0.15, 0.18, 0.21, 0.24),
    method = "t.test",
    label = "p.signif"
  )

# save
ggsave(
  glue(
    "{result_folder}Insulation_scores_200kb_windows_20kb_HiC_res-trimmed.pdf"
  ),
  plot = last_plot(),
  width = 8,
  height = 6
)

ggsave(
  glue(
    "{result_folder}Insulation_scores_200kb_windows_20kb_HiC_res-trimmed.png"
  ),
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 500
)

# TMPyP4, non-treated insulation score bigwigs (coming from 20 kb Hi-C)
tmpyp4 = "../data/Hi-C/TMPyP4_fastq_merge_insulation.bw"
nt = "../data/Hi-C/NT_fastq_merge_insulation.bw"

# TMPyP4 vs. non-treated IS scatter plot
scatter1 = plot_bw_bins_scatter(
  tmpyp4,
  nt, 
  bin_size = 2000,
  remove_top = 0.25,
  verbose = FALSE
)
scatter1
scatter1 = scatter1 + geom_point(size = 2, color = "#3182bd") +
  labs(title = "TMPyP4 vs. non-treated",
       x = "TMPyP4 aggr. IS",
       y = "non-treated aggr. IS") +
  guides(fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, color = "black",),
    axis.text.y = element_text(size = 10, color = "black")
  )
scatter1

# save
ggsave(
  glue(
    "{result_folder}Insulation_scores_NT-TMPyP4_scatter.pdf"
  ),
  plot = last_plot(),
  width = 8,
  height = 8
)

ggsave(
  glue(
    "{result_folder}Insulation_scores_NT-TMPyP4_scatter.png"
  ),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 500
)

tmpyp4_vs_nt = bw_bins(c(tmpyp4, nt),
                  bin_size = 2000,
                  remove_top = 0.1,
                  genome = "mm10")
tmpyp4_vs_nt = as_tibble(tmpyp4_vs_nt)

neg_IS = tmpyp4_vs_nt %>% dplyr::filter(TMPyP4_fastq_merge_insulation < 0 & NT_fastq_merge_insulation < 0) %>% 
  mutate(fc = TMPyP4_fastq_merge_insulation / NT_fastq_merge_insulation) %>% mutate(log_fc = log2(fc)) %>% 
  mutate(diff = abs(TMPyP4_fastq_merge_insulation - NT_fastq_merge_insulation))
high_diffs = neg_IS %>% dplyr::filter(diff >= 2) %>% 
  mutate(direction = ifelse(TMPyP4_fastq_merge_insulation < NT_fastq_merge_insulation, "DOWN", "UP"))

# negative IS changes 
bar = high_diffs %>% group_by(direction) %>% count() %>% ggplot(., aes(x = direction, y = n)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", fill = "#3182bd") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  labs(title = "negative IS changes, 2 kb bins - TMPyP4 vs. non-treated",
       x = "",
       y = "# of changes") +
  guides(fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 13, color = "black",),
    axis.text.y = element_text(size = 13, color = "black")
  )
bar

# save
ggsave(
  glue(
    "{result_folder}Neg_IS_changes-NT-TMPyP4.pdf"
  ),
  plot = last_plot(),
  width = 5,
  height = 5
)

ggsave(
  glue(
    "{result_folder}Neg_IS_changes-NT-TMPyP4.png"
  ),
  plot = last_plot(),
  width = 5,
  height = 5,
  dpi = 500
)


write_tsv(neg_IS, glue(
  "{result_folder}Insulation_scores-TMPyP4_vs_NT_negative_IS.tsv"
))




