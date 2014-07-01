#! /usr/bin/python

import os
import sys
import datetime

if len (sys.argv) < 2 or sys.argv[1] == '--help' or sys.argv[1] == '-h':
    sys.exit ('Usage: driveCem [ --dryrun | --commit | --ghost ]')

# See whether the friendly user really wants to include this model.
option = sys.argv[1]

# Open up some relevant files.
log = open ( './dat/log.txt', 'a' )
drv = open ( './dat/driver.txt', 'rw' )

# Print the date and time.
date = datetime.datetime.now()

# Read the driver.
params = []
for line in drv:
    if not ( line.startswith ('##') or line == '\n' ):
        
        fields = line.split ('= ')
        params.append ( fields[1].rstrip() )

# Tell the friendly user what he/she is about to get involved with.

print "\n\n --- You're set up for a(n) " + params[8] + " operation with " + params[4] + " degrees of rotation. ---\n\n"

print "Here's a summary of what you're trying to do: \n"
print "The directory you're going to grab mesh files from is: " + params[0]
print "The directory where your model files are from is: " + params[1]
print "The model type you've specified is: " + params[2]
print "The system of anisotropy is: " + params[3]
print "The angle of your requested rotation is (in degrees): " + params[4]
print "The rotation vector is: " + params[5] + " " + params[6] + " " + params[7]
print "You're attempting a " + params[8] + " operation."
print "The output symmetry system is: " + params[9]
print "Are you trying to overwrite the crust (interpolation only)?: " + params[10]
print "Are we treating what's hapenning as a kernel update? (interpolation only)?: " + params[11]

drv.close

# Warn the friendly user that it's possible they're opening up a bag of worms.
if option == '--commit':
    choice = raw_input ( 'Are you sure you want to go ahead? This will commit the parameters to the actual log file [y/n]: ' )
    if ( choice == 'n' ):
        sys.exit ( "K, I'm out." )
        
# If we want to commit, write the parameters to the log file.
if option == '--commit':

    drv = open ( './dat/driver.txt', 'r' )
    log.write ( '\n ----- ADDITION ----- \n' )
    log.write ( str (date) + '\n' )
    for line in drv:
        if not ( line.startswith ('##') or line == '\n' ):
            log.write (line)

elif option == '--ghost':

    print "I'm not keeping a record of this action, but I am going to make a sbatch script."

# Close files.
log.close

# Figure out what type of thing we're trying to do.
if ( params[8] == 'EXTRACT' and params[2] == 'SES3D' ):
    intent = 'extract_s3d'
if ( params[8] == 'EXTRACT' and params[2] == 'SPECFEM' ):
    intent = 'extract_spec'
if ( params[8] == 'CRUST' ):
    intent = 'add_crust'
if ( params[8] == 'INTERPOLATE' ):
    intent = 'interpolate'
if ( params[8] == 'TOPOGRAPHY' ):
    intent = 'topo'


# Write sbatch submission.
if option == '--commit' or option == '--ghost':

    slurm = open ( './jobs/job_' + intent + '.sbatch', 'w' )
    slurm.write ( '#!/bin/bash -l\n\n' )
    slurm.write ( '#SBATCH --nodes=1\n' )
    slurm.write ( '#SBATCH --ntasks=1\n' )
    slurm.write ( '#SBATCH --time=00:30:00\n' )
    slurm.write ( '#SBATCH --cpus-per-task=20\n' )
    slurm.write ( '#SBATCH --mem=32GB\n' )
    slurm.write ( '#SBATCH --exclusive\n' )
    slurm.write ( '#SBATCH --partition=fichtner_compute\n' )
    slurm.write ( '#SBATCH --output=jobs/job.' + intent + '.out\n' )
    slurm.write ( '#SBATCH --error=jobs/job.' + intent + '.err\n\n' )

    slurm.write ( 'module load intel netcdf\n' )

    slurm.write ( 'export OMP_NUM_THREADS=20\n' )
    slurm.write ( 'srun ' + os.environ['CEMHOME'] + '/bin/' + intent + '\n' )
    slurm.close ()

    print 'Created ./jobs/job_' + intent + '.sbatch. Submit that bad boy.\n'
