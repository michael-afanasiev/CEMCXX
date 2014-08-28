#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=32GB
#SBATCH --exclusive
#SBATCH --partition=fichtner_compute
#SBATCH --output=jobs/output/job.interpolate.%j.out
#SBATCH --error=jobs/output/job.interpolate.%j.err

module load intel netcdf
export OMP_NUM_THREADS=20
srun /users/afanasm/src/codeDevelopment/CEMCXX/bin/interpolate
