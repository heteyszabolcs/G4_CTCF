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

computeMatrix scale-regions -o ../results/deeptools/matrix_expr_lowmedium.mat.gz \
 -S ../data/CutNTag/bw/CTCF_AID_0h_merge_5million.norm.bw \
 ../data/CutNTag/bw/CTCF_AID_6h_merge.norm.bw \
 ../data/CutNTag/bw/G4_0h_AID_5million.norm.bw \
 ../data/CutNTag/bw/G4_6h_AID_5million.norm.bw \
 ../data/CutNTag/bw/G4_24h_AID_5million.norm.bw \
 -R ../results/rna_seq_deseq2/NT_lowmedium.bed \
 --regionBodyLength 2000 \
 --beforeRegionStartLength 3000 \
 --afterRegionStartLength 3000 \
 --skipZeros --missingDataAsZero \
 --samplesLabel "CTCF 0h AID" "CTCF 6h AID" "G4 0h AID" "G4 6h AID" "G4 24h AID" \

plotHeatmap -m ../results/deeptools/matrix_expr_lowmedium.mat.gz \
 -out "../results/deeptools/lowmediumexpr_scaleregion.pdf" \
 --heatmapHeight 14 \
 --colorMap "magma" \
 --yMin 0 \
 --yMax 10 \
 --zMax 10 \
 --yAxisLabel "" \
 --xAxisLabel "" \
 -z "low-medium expr. genes" \
 --outFileSortedRegions ../results/deeptools/matrix_expr_highmedium_sorted_regions.tsv \
 #--kmeans 2 \

plotHeatmap -m ../results/deeptools/matrix_expr_lowmedium.mat.gz \
 -out "../results/deeptools/lowmediumexpr_scaleregion.png" \
 --heatmapHeight 14 \
 --colorMap "magma" \
 --yMin 0 \
 --yMax 10 \
 --zMax 10 \
 --yAxisLabel "" \
 --xAxisLabel "" \
 -z "low-medium expr. genes" \
 #--kmeans 2 \