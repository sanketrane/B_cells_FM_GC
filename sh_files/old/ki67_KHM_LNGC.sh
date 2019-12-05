#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/FM_GC_ageCorrected
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G

srun Rscript scripts/ki67_KHM_LNGC.R
echo "Job Done!" | mail -A output/ki67_KHM_LNGC/job* sanketrn@gmail.com
