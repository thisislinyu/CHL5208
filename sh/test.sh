#!/bin/bash
#SBATCH --cpus-per-task 4
#SBATCH --time=0-00:30:00
#SBATCH --job-name test
#SBATCH --output=test_%j.txt
#SBATCH --mail-type=ALL

module load gcc/8.3.0 r/4.1.2

MAINDIR=/scratch/w/wxu/linyu/HDSI
cd ${MAINDIR}

Rscript ${MAINDIR}/R/05_top200_ccdb.R 

echo "Job complete!"
