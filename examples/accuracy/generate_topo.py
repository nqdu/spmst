import os
import numpy as np  
import matplotlib.pyplot as plt
import math

# define normalized 2D gaussian
def gaus2d_mountain(x=0, y=0, mx=0, my=0, sx=1, sy=1):
    return 1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))

def gaus2d_basin(x=0, y=0, mx=0, my=0, sx=1, sy=1):
    return (1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.))))*-1
### sx 和 sy 是 x 和 y 方向的传播，mx 和 my 是中心坐标。

n = 51
x = np.linspace(-5, 5,n)
y = np.linspace(-5, 5,n)
x, y = np.meshgrid(x, y) # get 2D variables instead of 1D


z = gaus2d_mountain(x, y,0,0,0.5,0.5)  ## mountain
# z = gaus2d_basin(x, y)  ## basin
# z=np.arctan(x)          ## plateau


## normal 0 ~ 5
zmax=np.max(z)
zmin=np.min(z)
for i in range(z.shape[1]):
    for j in range(z.shape[0]):
        z[i][j] = ((z[i][j]-zmin)/(zmax-zmin)) * 1000.

min_lon = x.min()
max_lon = x.max()
min_lat = y.min()
max_lat = y.max()
n_lons=x.shape[0]
n_lats=y.shape[0]
f = open('topo.txt','w')
f.write('%d %d\n' %(n_lons,n_lats))
f.write('%f %f %f %f\n'%(min_lon,max_lon,min_lat,max_lat))
for i in range(y.shape[0]):
    for j in range(x.shape[0]):
            f.write('%g\n' %(z[i,j])) 
f.close()