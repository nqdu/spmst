import numpy as np 
from scipy.interpolate import griddata

data = np.loadtxt("out.txt")
nlon = 265; nlat = 97
lon = np.linspace(89.2,90.3,nlon); lat = np.linspace(38.9,39.3,nlat)

x = np.zeros((nlon * nlat)); y = x.copy()

for i in range(nlat):
    for j in range(nlon):
        n = i * nlon + j
        x[n] = lon[j]; y[n] = lat[i]

z = griddata(data[:,0:2],data[:,2],(x,y),'nearest')

f = open("topo.dat","w")
f.write("265 97\n")
f.write("89.2 90.3 38.9 39.3\n")
for i in range(nlon*nlat):
    f.write("%f\n"%(z[i]))
f.close()