#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/FM_GC_ageCorrected
#SBATCH -N 4


srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/SPGC/ki67_TDT_SPGC.R "T1" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/SPGC/ki67_TDT_SPGC.R "T2" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/SPGC/ki67_TDD_SPGC.R "T1" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/SPGC/ki67_TDD_SPGC.R "T2" &
wait
