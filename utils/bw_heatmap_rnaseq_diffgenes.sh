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

computeMatrix scale-regions -o ../results/deeptools/matrix_diffgenes.mat.gz \
 -S ../data/CutNTag/bw/CTCF_AID_0h_merge_5million.norm.bw \
 ../data/CutNTag/bw/CTCF_AID_6h_merge.norm.bw \
 ../data/CutNTag/bw/G4_0h_AID_5million.norm.bw \
 ../data/CutNTag/bw/G4_6h_AID_5million.norm.bw \
 ../data/CutNTag/bw/G4_24h_AID_5million.norm.bw \
 -R ../results/rna_seq_deseq2/NT_low.bed \
 --regionBodyLength 2000 \
 --beforeRegionStartLength 3000 \
 --afterRegionStartLength 3000 \
 --skipZeros --missingDataAsZero \
 --samplesLabel "CTCF 0h AID" "CTCF 6h AID" "G4 0h AID" "G4 6h AID" "G4 24h AID" \

plotHeatmap -m ../results/deeptools/matrix_diffgenes.mat.gz \
 -out "../results/deeptools/diffgenes_fc0.25_p0.05.pdf" \
 --heatmapHeight 8 \
 --colorMap "magma" \
 --yMin 0 \
 --yMax 10 \
 --zMax 10 \
 -z "diff. expr. genes (log2FC: 0.25, adj.p: 0.05)" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_diffgenes.mat.gz \
 -out "../results/deeptools/diffgenes_fc0.25_p0.05.png" \
 --heatmapHeight 8 \
 --colorMap "magma" \
 --yMin 0 \
 --yMax 10 \
 --zMax 10 \
 -z "diff. expr. genes (log2FC: 0.25, adj.p: 0.05)" \
 --yAxisLabel "" \
 --xAxisLabel "" \
