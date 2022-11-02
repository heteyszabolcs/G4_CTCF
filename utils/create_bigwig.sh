#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 8:00:00
#SBATCH -M snowy
#SBATCH -J bigwigs

module load bioinfo-tools
module load R_packages
module load deepTools

Rscript create_bigwig.R
