#!/bin/bash -l
#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M rackham
#SBATCH -J deeptools

module load bioinfo-tools
module load deepTools

cd /proj/snic2020-6-3/SZABOLCS/HiC_G4s/utils

computeMatrix scale-regions -o ../results/deeptools/matrix_diffgenes_PQS.mat.gz \
 -S ../data/mm10_canPQS-regex_binary.bw \
 -R ../results/rna_seq_deseq2/RNA-Seq_treat_vs_contr_fc0.25_p0.05.bed \
 --regionBodyLength 2000 \
 --beforeRegionStartLength 3000 \
 --afterRegionStartLength 3000 \
 --skipZeros --missingDataAsZero \
 --samplesLabel "PQS" \

plotHeatmap -m ../results/deeptools/matrix_diffgenes_PQS.mat.gz \
 -out "../results/deeptools/diffgenes_fc0.25_p0.05_PQS.pdf" \
 --heatmapHeight 8 \
 --colorMap "binary" \
 -z "diff. expr. genes (log2FC: 0.25, adj.p: 0.05)" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_diffgenes_PQS.mat.gz \
 -out "../results/deeptools/diffgenes_fc0.25_p0.05_PQS.png" \
 --heatmapHeight 8 \
 --colorMap "binary" \
 -z "diff. expr. genes (log2FC: 0.25, adj.p: 0.05)" \
 --yAxisLabel "" \
 --xAxisLabel "" \
