#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=20
#SBATCH --partition=other_largemem
#SBATCH --output=jobs/job.interpolate.out
#SBATCH --error=jobs/job.interpolate.err

module load intel netcdf
export OMP_NUM_THREADS=20
srun /users/afanasm/src/codeDevelopment/CEMCXX/bin/interpolate
