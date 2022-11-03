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

computeMatrix reference-point -o ../results/deeptools/matrix.mat.gz \
 -S ../data/CutNTag/bw/CTCF_0h_AID_merge.bw \
 ../data/CutNTag/bw/CTCF_6h_AID_merge.bw \
 ../data/CutNTag/bw/G4_0h_AID_merge.bw \
 ../data/CutNTag/bw/G4_6h_AID_merge.bw \
 -R ../data/CutNTag/bed/CTCF_AID_0h_quant75_peaks.bed \
 -b 3000 -a 3000 --samplesLabel "CTCF 0h AID" "CTCF 6h AID" "G4 0h AID" "G4 6h AID"  \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix.mat.gz \
 -out "../results/deeptools/CTCF_G4_AID_treatment.pdf" \
 --refPointLabel "CTCF" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "CTCF peaks > 75th perc." \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix.mat.gz \
 -out "../results/deeptools/CTCF_G4_AID_treatment.png" \
 --refPointLabel "CTCF" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "CTCF peaks > 75th perc." \
 --yAxisLabel "" \
 --xAxisLabel "" \
