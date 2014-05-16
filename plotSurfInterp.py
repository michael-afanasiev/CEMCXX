#! /usr/local/bin/python

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

f = open ( './surfInterp.txt', 'r' )

x=[None]*5
y=[None]*5
d=[None]*5

lineNo = 1
for line in f:
  
  fields = line.split ()
  
  if lineNo == 1:
    x[0] = float (fields[0])
    y[0] = float (fields[1])
    d[0] = float (fields[2])
    
  if lineNo == 2:
    x[1] = float (fields[0])
    y[1] = float (fields[1])
    d[1] = float (fields[2])
  
  if lineNo == 3:
    x[2] = float (fields[0])
    y[2] = float (fields[1])
    d[2] = float (fields[2])
  
  if lineNo == 4:
    x[3] = float (fields[0])
    y[3] = float (fields[1])
    d[3] = float (fields[2])
  
  if lineNo == 5:
    x[4] = float (fields[0])
    y[4] = float (fields[1])
    d[4] = float (fields[2])
    
  lineNo = lineNo + 1
  
dMax = max(d)
dNorm = [i / dMax for i in d]
 
fig = plt.figure()
ax  = fig.add_subplot(111)

print dNorm
  
ax.scatter ( x, y, s=100, c=dNorm, marker='o', cmap=cm.hot )
plt.show()