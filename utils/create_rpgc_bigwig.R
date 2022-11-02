library("dplyr")
library("tidyverse")
library("data.table")

work_dir = "/proj/snic2020-6-3/SZABOLCS/HiC_G4s/data/ATAC-Seq/ATAC-Seq_mm10/"
bam_dir = "/proj/snic2020-6-3/SZABOLCS/HiC_G4s/data/ATAC-Seq/ATAC-Seq_mm10/results/bwa/mergedLibrary/"
bigwig_dir = "../results/rpgc_bigwig/"
# estimated size factors of Drosophila spike-ins
#size_f = fread("/proj/snic2020-6-3/SZABOLCS/HiC_G4s/data/ATAC-Seq/ATAC-Seq_mm10/spikein_size_factors.tsv")
bams = list.files(path = bam_dir, pattern = "bam$", full.names = TRUE)

for (bam in bams) {
  name = str_split(bam, ".mLb.clN.sorted.bam")[[1]][1]
  print(name)
  name = str_split(name, paste0(bam_dir, "/"))[[1]][2]
  print(name)
  # spikein_sf = size_f %>% filter(sample == name) %>% pull(size_factor) %>% as.numeric %>% round(., 2)
  print(name)
  system(paste0(
    "bamCoverage ",
    "-b ",
    bam,
    # " --scaleFactor ",
    # spikein_sf,
    " --normalizeUsing RPGC ",
    "--effectiveGenomeSize ",
    as.numeric(2652783500),
    " -o ",
    name,
    "_rpgc.bigwig"
  ))
  print(paste0(name, " is ready!"))
  
  # moving
  system(paste0("mv ", name, "_rpgc.bigwig ", bigwig_dir))
}
