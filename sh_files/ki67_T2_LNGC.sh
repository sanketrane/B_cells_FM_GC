#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/FM_GC_ageCorrected
#SBATCH -N 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr3506@cumc.columbia.edu

srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/ki67_SHM_LNGC.R "T2_data.csv" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/ki67_TDM_LNGC.R "T2_data.csv" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/ki67_SHM_SPGC.R "T2_data.csv" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/ki67_TDM_SPGC.R "T2_data.csv" &
wait
