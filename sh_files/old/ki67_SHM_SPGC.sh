#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/FM_GC_ageCorrected
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

srun Rscript scripts/ki67_SHM_SPGC.R
echo "Job Done!" | mail -A output/ki67_SHM_SPGC/job* sanketrn@gmail.com
