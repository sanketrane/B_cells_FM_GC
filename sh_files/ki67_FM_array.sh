#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%A.%a.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%A.%a.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/FM_GC_ageCorrected
#SBATCH --array=1-5
#SBATCH --cpus-per-task=25

srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/SPGC/ki67_SVM_SPGC.R "T2" 
