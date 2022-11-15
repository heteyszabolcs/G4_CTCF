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

computeMatrix reference-point -o ../results/deeptools/matrix_g4_kmeans_0h6h_AID.mat.gz \
 -S ../data/CutNTag/bw/CTCF_AID_0h_merge_5million.norm.bw \
 ../data/CutNTag/bw/CTCF_AID_6h_merge.norm.bw \
  ../data/CutNTag/bw/G4_0h_AID_5million.norm.bw \
 ../data/CutNTag/bw/G4_6h_AID_5million.norm.bw \
 ../data/CutNTag/bw/G4_24h_AID_5million.norm.bw \
 -R ../data/CutNTag/bed/G4_WT_peaks.narrowPeak \
 -b 3000 -a 3000 --samplesLabel "CTCF 0h AID" "CTCF 6h AID" "G4 0h AID" "G4 6h AID" "G4 24h AID" \
 --outFileNameMatrix ../results/deeptools/matrix_g4_kmeans_0h6h_AID.tsv \
 --outFileSortedRegions ../results/deeptools/matrix_g4_kmeans_0h6h_AID_sorted_regions.tsv \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_g4_kmeans_0h6h_AID.mat.gz \
 -out "../results/deeptools/G4_kmeans_0h_6hAID.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 --kmeans 2 \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --outFileSortedRegions ../results/deeptools/matrix_g4_kmeans_0h6h_AID_sorted_regions.tsv \

plotHeatmap -m ../results/deeptools/matrix_g4_kmeans_0h6h_AID.mat.gz \
 -out "../results/deeptools/G4_kmeans_0h_6hAID.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 --kmeans 2 \
 --yAxisLabel "" \
 --xAxisLabel "" \
