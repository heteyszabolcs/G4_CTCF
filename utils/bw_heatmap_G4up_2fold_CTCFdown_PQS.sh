#!/bin/bash -l
#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M snowy
#SBATCH -J deeptools

module load bioinfo-tools
module load deepTools

cd /proj/snic2020-6-3/SZABOLCS/HiC_G4s/utils

computeMatrix reference-point -o ../results/deeptools/matrix_G4up_2fold_CTCFdown_PQS.mat.gz \
 -S ../data/mm10_canPQS-regex_binary.bw \
 -R ../results/cutntag/G4up_2fold_CTCFdown_2fold_prom_regions.bed \
 --referencePoint center \
 --sortRegions keep \
 -b 10000 -a 10000 --samplesLabel "PQS"  \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_G4up_2fold_CTCFdown_PQS.mat.gz \
 -out "../results/deeptools/G4up_2fold_CTCFdown_upon6h_AID_PQS.pdf" \
 --sortRegions no \
 --refPointLabel "" \
 --heatmapHeight 14 \
 --colorMap "binary" \
 -z "G4 gained, CTCF loss promoters" \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --outFileSortedRegions ../results/deeptools/matrix_G4up_2fold_CTCFdown_PQS_sorted_regions.tsv \

plotHeatmap -m ../results/deeptools/matrix_G4up_2fold_CTCFdown_PQS.mat.gz \
 -out "../results/deeptools/G4up_2fold_CTCFdown_upon6h_AID_PQS.png" \
 --sortRegions no \
 --refPointLabel "" \
 --heatmapHeight 14 \
 --colorMap "binary" \
 -z "G4 gained, CTCF loss promoters" \
 --yAxisLabel "" \
 --xAxisLabel "" \
