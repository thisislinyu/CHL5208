#!/bin/bash
#SBATCH --nodes=20
#SBATCH --cpus-per-task 40
#SBATCH --time=0-12:00:00
#SBATCH --job-name sim_set2_20
#SBATCH --output=sim_set2_%j.txt
#SBATCH --mail-type=ALL

module load gcc/8.3.0 r/4.1.2

MAINDIR=/scratch/w/wxu/linyu/HDSI
cd ${MAINDIR}

Rscript ${MAINDIR}/R/99_sim_set2_ccdb20.R 

echo "Job complete!"
