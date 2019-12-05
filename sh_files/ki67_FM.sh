#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/FM_GC_ageCorrected
#SBATCH -N 2

srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/without_YFP_models/FM/ki67_TDT_FM.R "T1" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/without_YFP_models/FM/ki67_SHM_FM.R "T1" &
wait
