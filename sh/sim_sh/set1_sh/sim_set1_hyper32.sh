#!/bin/bash
#SBATCH --nodes=20
#SBATCH --cpus-per-task 40
#SBATCH --time=0-12:00:00
#SBATCH --job-name sim_set1_32
#SBATCH --output=sim_set1_%j.txt
#SBATCH --mail-type=ALL

module load gcc/8.3.0 r/4.1.2

MAINDIR=/scratch/w/wxu/linyu/HDSI
cd ${MAINDIR}

Rscript ${MAINDIR}/R/99_sim_set1_ccdb32.R 

echo "Job complete!"
