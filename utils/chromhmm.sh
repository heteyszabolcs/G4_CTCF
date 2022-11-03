#!/bin/bash -l
#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 16:00:00
#SBATCH -M snowy
#SBATCH -J wigglescout

module load bioinfo-tools
module load R_packages

# working directory
cd /proj/snic2020-6-3/SZABOLCS/HiC_G4s/utils/

Rscript /proj/snic2020-6-3/SZABOLCS/HiC_G4s/utils/chromhmm.R
