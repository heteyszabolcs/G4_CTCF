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

# computeMatrix reference-point -o ../results/deeptools/matrix_ctcf_kmeans.mat.gz \
 # -S ../data/CutNTag/bw/G4_0h_AID_5million.norm.bw \
 # ../data/CutNTag/bw/G4_6h_AID_5million.norm.bw \
 # ../data/CutNTag/bw/G4_24h_AID_5million.norm.bw \
 # ../data/CutNTag/bw/G4_48h_AID_5million.norm.bw \
 # ../data/CutNTag/bw/CTCF_AID_0h_merge_5million.norm.bw \
 # ../data/CutNTag/bw/CTCF_AID_6h_merge.norm.bw \
 # -R ../data/CutNTag/bed/CTCF_WT_peaks.narrowPeak \
 # -b 3000 -a 3000 --samplesLabel "G4 0h AID" "G4 6h AID" "G4 24h AID" "G4 48h AID" "CTCF 0h AID" "CTCF 6h AID" \
 # --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_ctcf_kmeans.mat.gz \
 -out "../results/deeptools/CTCF_kmeans.pdf" \
 --refPointLabel "CTCF" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 50 \
 --zMax 50 \
 --kmeans 3 \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --outFileSortedRegions ../results/deeptools/matrix_ctcf_kmeans_sorted_regions.tsv \

plotHeatmap -m ../results/deeptools/matrix_ctcf_kmeans.mat.gz \
 -out "../results/deeptools/CTCF_kmeans.png" \
 --refPointLabel "CTCF" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 50 \
 --zMax 50 \
 --kmeans 3 \
 --yAxisLabel "" \
 --xAxisLabel "" \
