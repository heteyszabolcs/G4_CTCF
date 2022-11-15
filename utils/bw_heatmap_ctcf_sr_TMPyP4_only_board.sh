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

computeMatrix scale-regions -o ../results/deeptools/matrix_ctcf_sr_NT_only.mat.gz \
 -S ../data/CutNTag/bw/CTCF_AID_0h_merge_5million.norm.bw \
 ../data/CutNTag/bw/CTCF_AID_6h_merge.norm.bw \
 -R ../data/Hi-C/bed/NT_only_fastq_merge_200kb_noheader.bed \
 --regionBodyLength 200000 \
 --beforeRegionStartLength 50000 \
 --afterRegionStartLength 50000 \
 --skipZeros --missingDataAsZero 
 #--samplesLabel "CTCF 0h AID" "CTCF 6h AID" \
 

plotHeatmap -m ../results/deeptools/matrix_ctcf_sr_NT_only.mat.gz \
 -out "../results/deeptools/CTCF_scaleregion_NT_only_board.pdf" \
 --startLabel "TAD boarder" \
 --endLabel "TAD boarder" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "TMPyP4 only boundaries" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_ctcf_sr_NT_only.mat.gz \
 -out "../results/deeptools/CTCF_scaleregion_NT_only_board.png" \
 --startLabel "TAD boarder" \
 --endLabel "TAD boarder" \
 --heatmapHeight 14 \
 --colorMap "YlOrRd" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "TMPyP4 only boundaries" \
 --yAxisLabel "" \
 --xAxisLabel "" \