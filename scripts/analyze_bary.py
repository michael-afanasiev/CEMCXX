#! /usr/local/bin/python

import matplotlib

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fName = 'Bary.txt'
f     = open (fName, 'r')

x     = [] 
y     = [] 
z     = []
point = []
l     = 0
first = True

# Set initial plot parameters.
fig = plt.figure()
ax  = fig.add_subplot (111, projection='3d')

ax.set_xlabel ('X-axis')
ax.set_ylabel ('Y-axis')
ax.set_zlabel ('Z-axis')
fig.suptitle  ("Testing for point inside Tet")

colours = ['green', 'red', 'blue', 'yellow', 'black', 'orange', 'purple']

plt.ion()
plt.show()

i = 0
for line in f:

  fields = line.split()
  if ( first ):
    point.append (float(fields[2]))
    point.append (float(fields[3]))
    point.append (float(fields[4]))
    first = False
    continue
    
  x.append (float(fields[2]))
  y.append (float(fields[3]))
  z.append (float(fields[4]))
  
  l = l+1
  
  if ( (l % 4) == 0 ):    
    
    first = True        

    xVert     = [x[0], x[1], x[2], x[3]]
    yVert     = [y[0], y[1], y[2], y[3]]
    zVert     = [z[0], z[1], z[2], z[3]]
    vertTuple = zip (xVert, yVert, zVert)

    vertConnec = [[1, 2, 0], [0, 2, 3], [0, 1, 3], [1, 2, 3]]
    poly3d     = [[vertTuple[vertConnec[ix][iy]] for iy in 
      range(len(vertTuple[0]))] for ix in range(len(vertTuple))]              

    
    polygon = Poly3DCollection (poly3d, facecolor=colours[i], linewidths=0.1,
      alpha=0.08)

    # Reset colorbar.
    i = i + 1
    if ( i == 7 ):
      i = 0
          
    ax.add_collection3d(polygon)
    ax.scatter (point[0], point[1], point[2], s=40, c='r')
        
    ax.set_xlim (point[0]-200, point[0]+200)
    ax.set_ylim (point[1]-200, point[1]+200)    
    ax.set_zlim (point[2]-200, point[2]+200)
    
    plt.draw()
    
    x = []
    y = []
    z = []
    
    check = raw_input()

f.close()
