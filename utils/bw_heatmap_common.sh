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

computeMatrix reference-point -o ../results/deeptools/matrix_comm.mat.gz \
 -S ../data/CutNTag/bw/CTCF_0h_AID_merge.bw \
 ../data/CutNTag/bw/CTCF_6h_AID_merge.bw \
 ../data/CutNTag/bw/CTCF_0h_TMPyP4_merge.bw \
 ../data/CutNTag/bw/CTCF_6h_TMPyP4_merge.bw \
 -R ../data/CutNTag/bed/common_CTCF_G4_AID_0h_quant75_peaks.bed \
 -b 3000 -a 3000 --samplesLabel "CTCF 0h AID" "CTCF 6h AID" "CTCF 0h TMPyP4" "CTCF 6h TMPyP4"  \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_comm.mat.gz \
 -out "../results/deeptools/common_CTCF_G4_both_treats.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "CTCF & G4" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_comm.mat.gz \
 -out "../results/deeptools/common_CTCF_G4_both_treats.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "CTCF & G4" \
 --yAxisLabel "" \
 --xAxisLabel "" \
