#!/bin/bash
#SBATCH --cpus-per-task 12
#SBATCH --time=0-12:00:00
#SBATCH --job-name sim_hyper_tune
#SBATCH --output=sim_hyper_tune_%j.txt
#SBATCH --mail-type=ALL

module load gcc/8.3.0 r/4.1.2

MAINDIR=/scratch/w/wxu/linyu/HDSI
cd ${MAINDIR}

Rscript ${MAINDIR}/R/06_simu_hyper_tune_ccdb.R 

echo "Job complete!"
