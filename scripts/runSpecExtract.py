#! /apps/monch/python-2.7.5/bin/python2.7

import collections
import os
import sys
import subprocess
import time

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

if ( len (sys.argv) < 5 or sys.argv[1] == '--help' ):
  print "Usage: \nNum Proc \n1st Proc \n1stReg \nlastReg"
  sys.exit()

numProc = int (sys.argv[1])
num1st  = int (sys.argv[2])
begReg  = int (sys.argv[3])
endReg  = int (sys.argv[4])

pwd = os.getcwd()

job_dict = collections.defaultdict(list)

cur_proc = 0
for i in range (num1st, numProc):
  job_dict[cur_proc].append(i)
  cur_proc += 1
  cur_proc %= size

for i in job_dict[rank]:
  for r in range (begReg-1, endReg):

    subprocess.call (pwd + '/bin/extract_spec ' + str (r) + ' ' + str (i),
    shell=True, stdout=None) 
