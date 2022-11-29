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

computeMatrix reference-point -o ../results/deeptools/matrix_G4up_2fold_CTCFdown.mat.gz \
 -S ../data/CutNTag/bw/G4_6h_AID_5million.norm.bw \
 ../data/CutNTag/bw/CTCF_AID_6h_merge.norm.bw \
 ../data/CutNTag/bw/G4_0h_AID_5million.norm.bw \
 ../data/CutNTag/bw/CTCF_AID_0h_merge_5million.norm.bw \
 -R ../results/cutntag/G4up_2fold_CTCFdown_2fold_prom_regions.bed \
 --referencePoint center \
 --sortRegions keep \
 -b 3000 -a 3000 --samplesLabel "G4 6h AID" "CTCF 6h AID" "G4 0h AID" "CTCF 0h AID"  \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_G4up_2fold_CTCFdown.mat.gz \
 -out "../results/deeptools/G4up_2fold_CTCFdown_upon6h_AID.pdf" \
 --sortRegions no \
 --refPointLabel "" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 15 \
 --zMax 15 \
 -z "G4 gained, CTCF loss promoters" \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --outFileSortedRegions ../results/deeptools/matrix_G4up_2fold_CTCFdown_sorted_regions.tsv \

plotHeatmap -m ../results/deeptools/matrix_G4up_2fold_CTCFdown.mat.gz \
 -out "../results/deeptools/G4up_2fold_CTCFdown_upon6h_AID.png" \
 --sortRegions no \
 --refPointLabel "" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 15 \
 --zMax 15 \
 -z "G4 gained, CTCF loss promoters" \
 --yAxisLabel "" \
 --xAxisLabel "" \
