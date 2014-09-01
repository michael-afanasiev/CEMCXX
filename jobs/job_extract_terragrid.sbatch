#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=10GB
#SBATCH --partition=fichtner_compute
#SBATCH --output=jobs/output/job.extract_spec.%j.out
#SBATCH --error=jobs/output/job.extract_spec.%j.err

module load intel netcdf
export OMP_NUM_THREADS=20
srun /users/afanasm/src/codeDevelopment/CEMCXX/bin/extract_terraGrid $1
