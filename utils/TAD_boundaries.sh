#!/bin/bash -l
#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 20:00:00
#SBATCH -M snowy
#SBATCH -J TAD_annot

# working directory
cd /proj/snic2020-6-3/SZABOLCS/HiC_G4s/utils/

module load bioinfo-tools
module load R_packages

Rscript /proj/snic2020-6-3/SZABOLCS/HiC_G4s/utils/TAD_boundaries.R