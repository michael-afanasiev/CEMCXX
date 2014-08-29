#!/bin/bash -l

#SBATCH --nodes=15
#SBATCH --ntasks=150
#SBATCH --ntasks-per-node=10
#SBATCH --time=06:00:00
#SBATCH --mem=32GB
#SBATCH --partition=fichtner_compute
#SBATCH --output=jobs/job.extract_spec.%j.out
#SBATCH --error=jobs/job.extract_spec.%j.err

module purge
module load mvapich2/1.9-intel python/2.7.5-intel-14.0.1 intel netcdf gcc
mpirun -n 150 /users/afanasm/src/codeDevelopment/CEMCXX/scripts/runSpecExtract.py 150 0 1 3
