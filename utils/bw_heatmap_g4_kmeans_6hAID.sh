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

computeMatrix reference-point -o ../results/deeptools/matrix_g4_kmeans_6hAID.mat.gz \
 -S ../data/CutNTag/bw/CTCF_6h_AID_merge.bw \
 -R ../data/CutNTag/bed/G4_WT_peaks.bed \
 --referencePoint center \
 -b 3000 -a 3000 --samplesLabel "CTCF 6h AID"  \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_g4_kmeans_6hAID.mat.gz \
 -out "../results/deeptools/G4_kmeans_6H_AID.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 --kmeans 2 \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_g4_kmeans_6hAID.mat.gz \
 -out "../results/deeptools/G4_kmeans_6H_AID.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 --kmeans 2 \
 --yAxisLabel "" \
 --xAxisLabel "" \
