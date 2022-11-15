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
# computeMatrix reference-point -o ../results/deeptools/matrix_ctcf_kmeans_24hAID.mat.gz \
 # -S ../data/CutNTag/bw/G4_24h_AID_merge.bw \
 # -R ../data/CutNTag/bed/CTCF_AID_0h_2_norm.bed \
 # -b 3000 -a 3000 --samplesLabel "G4 24 AID"  \
 # --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_ctcf_kmeans_24hAID.mat.gz \
 -out "../results/deeptools/CTCF_kmeans_24H_AID.pdf" \
 --refPointLabel "CTCF" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 200 \
 --zMax 200 \
 --kmeans 3 \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_ctcf_kmeans_24hAID.mat.gz \
 -out "../results/deeptools/CTCF_kmeans_24H_AID.png" \
 --refPointLabel "CTCF" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 200 \
 --zMax 200 \
 --kmeans 3 \
 --yAxisLabel "" \
 --xAxisLabel "" \
