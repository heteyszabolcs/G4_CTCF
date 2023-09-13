if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "wigglescout",
               "ggpubr",
               "glue",
               "RColorBrewer",
               "rstatix",
               "effsize")

# output folder
result_folder = "../results/wigglescout/"

# Minute-ChIP bigwigs
bws = list.files("../data/Minute-ChIP/",
                 pattern = "*scaled.bw",
                 full.names = TRUE)
bws = bws[grep("unscaled", bws, invert = TRUE)]
bws = bws[grep("pooled", bws, invert = FALSE)]

## helper functions
# create scatter plot between two Minute-ChIP bigwigs
create_scatter = function(bigwig1, bigwig2) {
  name1 = strsplit(strsplit(bigwig1, "../data/Minute-ChIP/")[[1]][2],
                   "_pooled.mm10.scaled.bw")[[1]][1]
  name2 = strsplit(strsplit(bigwig2, "../data/Minute-ChIP/")[[1]][2],
                   "_pooled.mm10.scaled.bw")[[1]][1]
  
  p = plot_bw_loci_scatter(
    bigwig1,
    bigwig2,
    loci =
      "../data/CutNTag_ChIP-Seq/bed/merged-CTCFonly_CTCFG4_common_with_promoters.bed",
    highlight = "../data/CutNTag_ChIP-Seq/bed/CTCF_G4_common_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed",
    highlight_colors = "#de2d26",
    remove_top = 0.001,
    verbose = FALSE
  ) + 
    xlim(0, 40) +
    ylim(0, 40) +
    theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black")
  ) +
    guides(color = "none", fill = "none") +
    geom_abline(slope = 1,
                intercept = 0,
                color = "#de2d26",
                linetype = "dashed") +
    labs(title = "G4-CTCF (red) & CTCF only (grey) promoters",
         x = name1,
         y = name2)
  print(p)
  # save
  ggsave(
    glue("{result_folder}{name1}-{name2}_G4CTCFproms-CTCFG4_red.pdf"),
    plot = p,
    width = 3,
    height = 3
  )
  return(p)
}

# return mean aggr. CTCF signals over a bed of interest
aggr_signal = function(bigwigs, bed, compound, regions) {
  
  output = bw_loci(
    bigwigs,
    loci = bed
  )
  output = as.data.frame(output)
  output = output %>% select(-seqnames,-start,-end,-width,-strand)
  colnames(output) = c("non-treated", paste0("6h ", compound), paste0("24h ", compound))
  output = output %>% pivot_longer("non-treated":paste0("24h ", compound),
                                                   names_to = "treatment",
                                                   values_to = "CTCF_signal") %>% 
    mutate(regions = regions)
  
  return(output)
}

# return boxplot about the mean aggr. data
make_aggr_boxplot = function(bp_input, compound = "TMPyP4", prom_stat = "prom") {
  
  # statistical tests
  stat.test = bp_input %>%
    group_by(treatment) %>%
    t_test(CTCF_signal ~ regions) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>% 
    mutate(stat = "t-test")
  print(stat.test)
  stat.test2 = bp_input %>%
    group_by(regions) %>%
    t_test(CTCF_signal ~ treatment) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>% 
    mutate(stat = "t-test")
  print(stat.test2)
  
  stat.test_wilcoxon = bp_input %>%
    group_by(treatment) %>%
    rstatix::wilcox_test(CTCF_signal ~ regions) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>% 
    mutate(stat = "wilcoxon")
  print(stat.test_wilcoxon)
  
  # Cohen's d effect sizes
  cohen_d_regions = cohen.d(d = bp_input$CTCF_signal, f = bp_input$regions)
  print(glue(
    "level of regions, {cohen_d_regions$method}: {as.character(round(cohen_d_regions$estimate, 2))}"
  ))
  
  trt24 = paste0("24h ", compound)
  cohen_input = bp_input %>% dplyr::filter(treatment == trt24 |
                                             treatment == "non-treated")
  cohen_d_24trt = cohen.d(d = cohen_input$CTCF_signal, f = cohen_input$treatment)
  print(glue(
    "level of 24h {compound} treatment, {cohen_d_24trt$method}: {as.character(round(cohen_d_24trt$estimate, 2))}"
  ))
  trt6 = paste0("6h ", compound)
  cohen_input = bp_input %>% dplyr::filter(treatment == trt6 |
                                             treatment == "non-treated")
  cohen_d_6trt = cohen.d(d = cohen_input$CTCF_signal, f = cohen_input$treatment)
  print(glue(
    "level of 6h {compound} treatment, {cohen_d_6trt$method}: {as.character(round(cohen_d_6trt$estimate, 2))}"
  ))
  
  # plot
  p = ggplot(bp_input, aes(x = regions, y = CTCF_signal)) +
    geom_boxplot(alpha = 1, aes(fill = treatment)) +
    scale_fill_manual(values = c("#de2d26", "#fc9272", "#bdbdbd")) +
    ylim(0, 150) +
    labs(
      title = glue("aggr. CTCF signal over {prom_stat}"),
      x = "",
      y = "mean aggr. CTCF signal",
      fill = ""
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(size = 13),
      axis.text.x = element_text(size = 15, color = "black", angle = 0),
      axis.title.y = element_text(size = 15, color = "black"),
      axis.text.y = element_text(size = 15, color = "black")
    )
  p
  
  stat.test = stat.test %>% add_xy_position(x = "regions", dodge = 0.75, group = "treatment")
  p = p + stat_pvalue_manual(
    stat.test, y.position = 135, tip.length = 0.01,
    bracket.nudge.y = -12, step.increase = c(0.08, 0.08, 0.08),
    label = "p.adj.signif"
  ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(0, 150))
  
  stat.test2 = stat.test2 %>% add_xy_position(x = "regions", dodge = 0.78, 
                                              group = "treatment")
  p = p + stat_pvalue_manual(
    stat.test2, y.position = 85, tip.length = 0.01,
    bracket.nudge.y = -12, step.increase = c(0.06, 0.06, 0.06),
    label = "p.adj.signif"
  ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(0, 150))
  print(p)
  
  ggsave(
    glue("{result_folder}aggr_CTCF_{prom_stat}_signal-{compound}_trt.pdf"),
    plot = p,
    width = 6,
    height = 4
  )
  
  return(p)
  
}

# execution
# loop through all pooled & scaled bigwig file
# lapply(bws, function(x) {
#   lapply(bws, function(y) {
#     create_scatter(x, y)
#   })
# })


### with promoter bed
## TMPyP4 treatment - CTCF Minute-ChIP
tmpyp4_CTCF_G4 = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_TMPyP4_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_TMPyP4_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_G4_common_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "TMPyP4", regions = "CTCF-G4")

tmpyp4_CTCFonly = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_TMPyP4_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_TMPyP4_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_no_G4_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "TMPyP4", regions = "CTCF only")

tmpyp4_bp_input = rbind(tmpyp4_CTCFonly, tmpyp4_CTCF_G4)

tmpyp4_bp = make_aggr_boxplot(bp_input = tmpyp4_bp_input, compound = "TMPyP4")

## PDS treatment - CTCF Minute-ChIP
pds_CTCF_G4 = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PDS_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PDS_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_G4_common_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "PDS", regions = "CTCF-G4")

pds_CTCFonly = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PDS_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PDS_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_no_G4_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "PDS", regions = "CTCF only")

pds_bp_input = rbind(pds_CTCFonly, pds_CTCF_G4)

pds_bp = make_aggr_boxplot(bp_input = pds_bp_input, compound = "PDS")

## PhenDC treatment - CTCF Minute-ChIP
phen_CTCF_G4 = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PhenDC3_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PhenDC3_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_G4_common_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "PhenD3", regions = "CTCF-G4")

phen_CTCFonly = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PhenDC3_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PhenDC3_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_no_G4_WT_with_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "PhenD3", regions = "CTCF only")

phen_bp_input = rbind(phen_CTCFonly, phen_CTCF_G4)

phen_bp = make_aggr_boxplot(bp_input = phen_bp_input, compound = "PhenD3")

# arrange
ggarrange(tmpyp4_bp, pds_bp, phen_bp, ncol = 2, nrow = 2)
ggsave(
  glue("{result_folder}aggr_CTCF_prom_signal-TMPyP4_PDS_PhenDC3_trt.pdf"),
  plot = last_plot(),
  width = 12,
  height = 8
)

### without promoter bed
## TMPyP4 treatment - CTCF Minute-ChIP
tmpyp4_CTCF_G4 = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_TMPyP4_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_TMPyP4_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_G4_common_WT_without_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "TMPyP4", regions = "CTCF-G4")

tmpyp4_CTCFonly = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_TMPyP4_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_TMPyP4_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_no_G4_WT_without_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "TMPyP4", regions = "CTCF only")

tmpyp4_bp_input = rbind(tmpyp4_CTCFonly, tmpyp4_CTCF_G4)

tmpyp4_bp = make_aggr_boxplot(bp_input = tmpyp4_bp_input, compound = "TMPyP4", prom_stat = "non_prom")

## PDS treatment - CTCF Minute-ChIP
pds_CTCF_G4 = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PDS_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PDS_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_G4_common_WT_without_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "PDS", regions = "CTCF-G4")

pds_CTCFonly = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PDS_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PDS_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_no_G4_WT_without_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "PDS", regions = "CTCF only")

pds_bp_input = rbind(pds_CTCFonly, pds_CTCF_G4)

pds_bp = make_aggr_boxplot(bp_input = pds_bp_input, compound = "PDS", prom_stat = "non_prom")

## PhenDC treatment - CTCF Minute-ChIP
phen_CTCF_G4 = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PhenDC3_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PhenDC3_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_G4_common_WT_without_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "PhenD3", regions = "CTCF-G4")

phen_CTCFonly = aggr_signal(c(
  "../data/Minute-ChIP/CTCF_E14_NT_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PhenDC3_6h_pooled.mm10.scaled.bw",
  "../data/Minute-ChIP/CTCF_E14_PhenDC3_24h_pooled.mm10.scaled.bw"
), bed = "../data/CutNTag_ChIP-Seq/bed/CTCF_no_G4_WT_without_promoters_missingaszero_skipzero.ordered.regions.bed",
compound = "PhenD3", regions = "CTCF only")

phen_bp_input = rbind(phen_CTCFonly, phen_CTCF_G4)

phen_bp = make_aggr_boxplot(bp_input = phen_bp_input, compound = "PhenD3", prom_stat = "non_prom")

# arrange
ggarrange(tmpyp4_bp, pds_bp, phen_bp, ncol = 2, nrow = 2)
ggsave(
  glue("{result_folder}aggr_CTCF_nonprom_signal-TMPyP4_PDS_PhenDC3_trt.pdf"),
  plot = last_plot(),
  width = 12,
  height = 8
)
