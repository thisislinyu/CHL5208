#!/bin/bash
#SBATCH --cpus-per-task 24
#SBATCH --time=0-12:00:00
#SBATCH --job-name real_hyper_tune
#SBATCH --output=real_hyper_tune_%j.txt
#SBATCH --mail-type=ALL

module load gcc/8.3.0 r/4.1.2

MAINDIR=/scratch/w/wxu/linyu/HDSI
cd ${MAINDIR}

Rscript ${MAINDIR}/R/07_real_dat_hyper_tune.R 

echo "Job complete!"
