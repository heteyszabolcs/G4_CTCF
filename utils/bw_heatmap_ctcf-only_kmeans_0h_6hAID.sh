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

computeMatrix reference-point -o ../results/deeptools/matrix_ctcf-only_kmeans_0h6h_AID.mat.gz \
 -S ../data/CutNTag/bw/CTCF_AID_0h_merge_5million.norm.bw \
 ../data/CutNTag/bw/CTCF_AID_6h_merge.norm.bw \
 ../data/CutNTag/bw/G4_0h_AID_5million.norm.bw \
 ../data/CutNTag/bw/G4_6h_AID_5million.norm.bw \
 ../data/CutNTag/bw/G4_24h_AID_5million.norm.bw \
 -R ../data/CutNTag/bed/CTCF_no_G4_WT_peaks.bed \
 --referencePoint center \
 -b 3000 -a 3000 --samplesLabel "CTCF 0h AID" "CTCF 6h AID" "G4 0h AID" "G4 6h AID" "G4 24h AID" \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_ctcf-only_kmeans_0h6h_AID.mat.gz \
 -out "../results/deeptools/CTCF-only_kmeans_0h_6hAID.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 20 \
 --zMax 20 \
 --kmeans 2 \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --outFileSortedRegions ../results/deeptools/matrix_ctcf-only_kmeans_0h6h_AID.tsv \

plotHeatmap -m ../results/deeptools/matrix_ctcf-only_kmeans_0h6h_AID.mat.gz \
 -out "../results/deeptools/CTCF-only_kmeans_0h_6hAID.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 20 \
 --zMax 20 \
 --kmeans 2 \
 --yAxisLabel "" \
 --xAxisLabel "" \
