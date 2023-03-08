# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("ggrepel")
  library("HiCcompare")
})


# Bonev et al. chr1 HiC data
chr1_norm = cooler2bedpe(path = "../data/Bonev_ES_10k_chr1_norm.cool")
head(chr1_norm$cis$chr1)

# Bonev et al. chr1 normalized (HiCExplorer, hicNormalize function - norm_range mode)
chr1 = cooler2bedpe(path = "../data/Bonev_ES_10k_chr1.cool")
head(chr1$cis$chr1)

# create input for meta-waffle
a = as_tibble(chr1$cis$chr1)
b = as_tibble(chr1_norm$cis$chr1)
a = a %>% unite(id, chr2, start2, end2, sep = "_", remove = FALSE) %>% distinct(., id, .keep_all = TRUE) 
b = b %>% unite(id, chr2, start2, end2, sep = "_", remove = FALSE) %>% distinct(., id, .keep_all = TRUE) 
c = a %>% inner_join(., b, by = "id")
d = tibble(start = c$start2.x, end = c$end2.x, raw = round(c$IF.x, 3), norm = round(c$IF.y, 3))

write_tsv(d, "../data/Bonev_et_al_chr1_norm.tsv", col_names = FALSE)
