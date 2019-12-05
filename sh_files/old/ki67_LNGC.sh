#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/FM_GC_ageCorrected
#SBATCH -N 3

srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/ki67_TDM_LNGC.R "FM_data.csv" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/ki67_TDM_LNGC.R "T2_data.csv" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/ki67_TDM_LNGC.R "T1_data.csv" &
wait
