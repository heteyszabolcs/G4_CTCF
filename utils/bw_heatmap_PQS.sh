#!/bin/bash -l
#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M snowy
#SBATCH -J deeptools

module load bioinfo-tools
module load deepTools

cd /proj/snic2020-6-3/SZABOLCS/HiC_G4s/utils

computeMatrix reference-point -o ../results/deeptools/matrix_G4_no_CTCF-PQS.mat.gz \
 -S ../data/bw/mm10_canPQS-regex_binary.bw \
 -R ../data/CutNTag/bed/G4_no_CTCF_WT_peaks.bed \
 --referencePoint center \
 -b 3000 -a 3000 --samplesLabel "PQS"  \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_G4_no_CTCF-PQS.mat.gz \
 -out "../results/deeptools/G4_no_CTCF-PQS.pdf" \
 --refPointLabel "" \
 --heatmapHeight 14 \
 --colorMap "binary" \
 -z "G4 only peaks" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_G4_no_CTCF-PQS.mat.gz \
 -out "../results/deeptools/G4_no_CTCF-PQS.png" \
 --refPointLabel "" \
 --heatmapHeight 14 \
 --colorMap "binary" \
 -z "G4 only peaks" \
 --yAxisLabel "" \
 --xAxisLabel "" \
