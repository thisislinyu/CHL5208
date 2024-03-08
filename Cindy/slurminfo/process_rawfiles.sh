#!/bin/bash

CHR=$1

MAINDIR=~/project/cim1/fracture_coxph/recoded_geno
cd ${MAINDIR}

module load r/3.4.0

Rscript ${MAINDIR}/process_rawfile_1.R $CHR

# NOTES
#This is an example R job submission script.
#You can run everything starting from line 5 straight into the "login terminal"!
#First line must be included to set the programming language.
#CHR is an "argument" that is passed along to your R script (like setting CHR=22, then every time you refer to CHR in the Rscript, that means CHR=22).
#MAINDIR is a variable you set in bash; this is the directory I want the R script to run from.
#You must load R (module load).
#Rscript is the slurm command to run an R script.
#Submit job command lines are below; the first line is needed to set the "argument".
#CHR=22
#sbatch -J results_${CHR} --mem=8G --time=00-02:00:00 --account=def-yyasui --wrap="${MAINDIR}/process_rawfile_1.sh ${sample}"

#Important: You must convert .sh to unix and make .sh scripts executable before you try to submit the .sh job. See below.
#dos2unix process_rawfiles.sh
#chmod 770 process_rawfiles.sh