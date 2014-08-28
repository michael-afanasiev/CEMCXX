#! /usr/bin/python

f = open ('./submitTerraGrid.sh', 'w')

for i in range (256):
  f.write ('sbatch ./jobs/job_extract_terragrid.sbatch ' + str(i) + '\n')
  
f.close()

