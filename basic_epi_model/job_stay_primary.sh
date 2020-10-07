#!/bin/bash

#SBATCH --job-name=sp_orig
#SBATCH --output=sp_orig.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

module load R/3.5.3-foss-2018a-X11-20180131
R --slave -f epidemia_stay_primary.R

