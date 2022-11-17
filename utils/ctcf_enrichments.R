suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("glue")
  library("wigglescout")
  library("GenomeInfoDb")
})

# annotation function
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# folders
result_folder = "../results/cutntag/"
bigwig_folder = "../data/CutNTag/bw/"

# for subsetting
peaks = "../data/CutNTag/bed/CTCF_AID_0h_quant75_peaks.bed"

calc_fold_changes = function(bigwig1,
                             bigwig2,
                             bed) {
  # fold changes by wigglescout
  fc = bw_loci(bigwig1,
               bigwig2,
               loci = bed)
  fc = as.data.frame(fc)
  fc_colname = colnames(fc)[6]
  
  fc = mm10_annotation(
    fc,
    seqname_col = "seqnames",
    start_col = "start",
    end_col = "end",
    feature_1 = fc_colname,
    feature_2 = fc_colname
  )
  fc = fc %>% dplyr::select(seqnames,
                            start,
                            end,
                            feature_1,
                            annotation,
                            distanceToTSS,
                            gene = SYMBOL)
  return(fc)
  
  
}

bigwigs = list.files(bigwig_folder, full.names = TRUE)


fc_annot = calc_fold_changes(bigwig1 = "../data/CutNTag/bw/CTCF_0h_AID_merge.bw",
                             bigwig2 = "../data/CutNTag/bw/CTCF_6h_AID_merge.bw",
                             bed = peaks)
fc_tmpyp4_annot = calc_fold_changes(bigwig1 = "../data/CutNTag/bw/CTCF_0h_TMPyP4_merge.bw",
                                    bigwig2 = "../data/CutNTag/bw/CTCF_6h_TMPyP4_merge.bw",
                                    bed = peaks)

# violin plot visualization of genomic element fold changes
vis = fc_annot %>%
  dplyr::rename(CTCF_0h_vs_CTCF_6h_AID = feature_1) %>% 
  dplyr::filter(!CTCF_0h_vs_CTCF_6h_AID == Inf) %>% 
  dplyr::filter(CTCF_0h_vs_CTCF_6h_AID < quantile(CTCF_0h_vs_CTCF_6h_AID, 0.9)) %>%
  mutate(annotation_v2 = ifelse(str_detect(annotation, "Exon"), "Exon", annotation)) %>%
  mutate(annotation_v2 = ifelse(str_detect(annotation_v2, "Intron"), "Intron", annotation_v2))  %>%
  ggplot(.,
         aes(x = annotation_v2, y = CTCF_0h_vs_CTCF_6h_AID, fill = annotation_v2)) +
  geom_violin(trim = TRUE, color = "black") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "CTCF 0h AID vs. CTCF 6h AID",
       x = " ",
       y = "fold change",
       fill = " ") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 13,
      color = "black",
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(size = 13, color = "black")
  )
vis

ggsave(
  glue("{result_folder}bw_loci_CTCF_AID_violins.pdf"),
  plot = vis,
  width = 8,
  height = 5,
  device = "pdf"
)

ggsave(
  glue("{result_folder}bw_loci_CTCF_AID_violins.png"),
  plot = vis,
  width = 8,
  height = 5,
  dpi = 300
)

tmpyp4_vis = fc_tmpyp4_annot %>%
  dplyr::rename(CTCF_0h_vs_CTCF_6h_TMPyP4 = feature_1) %>% 
  dplyr::filter(!CTCF_0h_vs_CTCF_6h_TMPyP4 == Inf) %>% 
  dplyr::filter(CTCF_0h_vs_CTCF_6h_TMPyP4 < 50) %>%
  mutate(annotation_v2 = ifelse(str_detect(annotation, "Exon"), "Exon", annotation)) %>%
  mutate(annotation_v2 = ifelse(str_detect(annotation_v2, "Intron"), "Intron", annotation_v2))  %>%
  ggplot(.,
         aes(x = annotation_v2, y = CTCF_0h_vs_CTCF_6h_TMPyP4, fill = annotation_v2)) +
  geom_violin(trim = TRUE, color = "black") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "CTCF 0h TMPyP4 vs. CTCF 6h TMPyP4",
       x = " ",
       y = "fold change",
       fill = " ") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 13,
      color = "black",
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(size = 13, color = "black")
  )
tmpyp4_vis

ggsave(
  glue("{result_folder}bw_loci_CTCF_TMPyP4_violins.pdf"),
  plot = tmpyp4_vis,
  width = 8,
  height = 5,
  device = "pdf"
)

ggsave(
  glue("{result_folder}bw_loci_CTCF_TMPyP4_violins.png"),
  plot = tmpyp4_vis,
  width = 8,
  height = 5,
  dpi = 300
)

# with the significant G4 peaks
peaks = "../data/CutNTag/bed/G4_AID_0h_quant75_peaks.bed"

fc_g4_aid_annot = calc_fold_changes(bigwig1 = "../data/CutNTag/bw/G4_0h_AID_merge.bw",
                             bigwig2 = "../data/CutNTag/bw/G4_48h_AID_merge.bw",
                             bed = peaks)

# violin plot visualization of genomic element fold changes
g4_aid_vis = fc_g4_aid_annot %>%
  dplyr::rename(G4_0h_vs_G4_48h_AID = feature_1) %>% 
  dplyr::filter(!G4_0h_vs_G4_48h_AID == Inf) %>% 
  #dplyr::filter(G4_0h_vs_G4_48h_AID < quantile(G4_0h_vs_G4_48h_AID, 0.9)) %>%
  mutate(annotation_v2 = ifelse(str_detect(annotation, "Exon"), "Exon", annotation)) %>%
  mutate(annotation_v2 = ifelse(str_detect(annotation_v2, "Intron"), "Intron", annotation_v2))  %>%
  ggplot(.,
         aes(x = annotation_v2, y = G4_0h_vs_G4_48h_AID, fill = annotation_v2)) +
  geom_violin(trim = FALSE, color = "black") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "G4 0h AID vs. G4 48h AID",
       x = " ",
       y = "fold change",
       fill = " ") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 13,
      color = "black",
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(size = 13, color = "black")
  )
g4_aid_vis

ggsave(
  glue("{result_folder}bw_loci_G4_AID_violins.pdf"),
  plot = g4_aid_vis,
  width = 8,
  height = 5,
  device = "pdf"
)

ggsave(
  glue("{result_folder}bw_loci_G4_AID_violins.png"),
  plot = g4_aid_vis,
  width = 8,
  height = 5,
  dpi = 300
)

