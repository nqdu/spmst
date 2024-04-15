import os
import numpy as np  
import matplotlib.pyplot as plt
import math

# define normalized 2D gaussian
def gaus2d_mountain(x=0, y=0, mx=0, my=0, sx=1, sy=1):
    return 1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))

def gaus2d_basin(x=0, y=0, mx=0, my=0, sx=1, sy=1):
    return - gaus2d_mountain(x,y,mx,my,sx,sy)
    #return (1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.))))*-1

# input args
nlon,nlat = 51,51
xmin, xmax = -5,5
ymin, ymax = -5,5
max_height = 1

# compute topography
x = np.linspace(xmin, xmax,nlon)
y = np.linspace(ymin, ymax,nlat)
X, Y = np.meshgrid(x, y) # get 2D variables instead of 1D
z = gaus2d_mountain(X,Y,0,0,0.5,0.5)  ## mountain

# normalize
z = z / np.max(np.abs(z)) * max_height

f = open('topo.txt','w')
f.write('%d %d\n' %(nlon,nlat))
f.write('%f %f\n%f %f\n'%(xmin,xmax,ymin,ymax))
for i in range(nlat):
    for j in range(nlon):
            f.write('%g\n' %(z[i,j])) 
f.close()