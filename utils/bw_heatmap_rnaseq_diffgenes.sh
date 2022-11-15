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

computeMatrix reference-point -o ../results/deeptools/matrix_diffgenes.mat.gz \
 -S ../data/CutNTag/bw/G4_0h_AID_merge.bw \
 ../data/CutNTag/bw/CTCF_0h_AID_merge.bw \
 -R ../results/rna_seq_deseq2/RNA-Seq_treat_vs_contr_fc0.25_p0.05.bed \
 -b 3000 -a 3000 --samplesLabel "G4 0h AID" "CTCF 0h AID"  \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_diffgenes.mat.gz \
 -out "../results/deeptools/diffgenes_fc0.25_p0.05.pdf" \
 --refPointLabel "TSS" \
 --heatmapHeight 8 \
 --colorMap "magma" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "diff. expr. genes (log2FC: 0.25, adj.p: 0.05)" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_diffgenes.mat.gz \
 -out "../results/deeptools/diffgenes_fc0.25_p0.05.png" \
 --refPointLabel "TSS" \
 --heatmapHeight 8 \
 --colorMap "magma" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "diff. expr. genes (log2FC: 0.25, adj.p: 0.05)" \
 --yAxisLabel "" \
 --xAxisLabel "" \
