#! /usr/local/bin/python

import os
import sys
import subprocess

if ( len (sys.argv) < 5 or sys.argv[1] == '--help' ):
  print "Usage: \nNum Proc \n1st Proc \n1stReg \nlastReg"
  sys.exit()

numProc = int (sys.argv[1])
num1st  = int (sys.argv[2])
begReg  = int (sys.argv[3])
endReg  = int (sys.argv[4])

pwd = os.getcwd()

# numProc = input ('Enter number of processors: ')
# num1st  = input ('Enter first processor: ')
# begReg, endReg = input ('Enter first, last region: ')


for r in range (begReg-1, endReg):
  for i in range (num1st, numProc):
    
    subprocess.call (pwd + '/bin/extract_spec ' + str (r) + ' ' + str (i),
    shell=True, stdout=None)
  
