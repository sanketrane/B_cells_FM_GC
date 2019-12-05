#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/FM_GC_ageCorrected
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH -N 4

srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript scripts/ki67_SHM_LNGC.R &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript scripts/ki67_KHM_LNGC.R &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript scripts/ki67_TDM_LNGC.R &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript scripts/ki67_INC_LNGC.R &
wait
