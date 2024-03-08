#!/bin/bash
#SBATCH --nodes=20
#SBATCH --cpus-per-task 40
#SBATCH --time=0-08:00:00
#SBATCH --job-name sim_set1_4test
#SBATCH --output=sim_set1_%j.txt
#SBATCH --mail-type=ALL

module load gcc/8.3.0 r/4.1.2

MAINDIR=/scratch/w/wxu/linyu/HDSI
cd ${MAINDIR}

Rscript ${MAINDIR}/R/99_sim_set1_ccdb4_test.R 

echo "Job complete!"