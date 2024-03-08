#!/bin/bash
#SBATCH --cpus-per-task 4
#SBATCH --time=0-1:00:00
#SBATCH --job-name testplot
#SBATCH --output=test_output_%j.txt
#SBATCH --mail-type=ALL

echo 'Hello, world!'
