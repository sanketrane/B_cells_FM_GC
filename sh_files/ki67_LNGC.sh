#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/FM_GC_ageCorrected
#SBATCH -N 6

srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/LNGC/ki67_TDT_LNGC.R "T1" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/LNGC/ki67_TDT_LNGC.R "T2" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/LNGC/ki67_TDT_LNGC.R "FM" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/LNGC/ki67_TDD_LNGC.R "T1" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/LNGC/ki67_TDD_LNGC.R "T2" &
srun -N 1 --ntasks=1 --ntasks-per-node=1 Rscript --vanilla scripts/YFP_informed_models/LNGC/ki67_TDD_LNGC.R "FM" &
wait
