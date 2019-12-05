#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --job-name=Test
#SBATCH --output=slurm_%a.out
/usr/lib/R/bin/Rscript --vanilla slurm_run.R
