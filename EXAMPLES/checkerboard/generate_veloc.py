import numpy as np
import matplotlib.pyplot as plt 

lonmin,latmin = 89.2,  38.9 # lonmin latmin
lonmax,latmax = 90.3,  39.3  # lonmax latmax
nlon = 40; nlat = 40
lon = np.linspace(lonmin,lonmax,nlon)
lat = np.linspace(latmin,latmax,nlat)

veloc = np.zeros((nlat,nlon))  + 3.0

f = open("velocinit.in","w")
for i in range(nlat):
    for j in range(nlon):
        f.write("%f\n"%veloc[i,j])
f.close()

x = np.arange(1,nlon+1) - 1
y = np.arange(1,nlat+1) - 1
veloc1 = veloc * 1.0
for i in range(nlat):
    for j in range(nlon):
        term = np.sin(x[j] * 0.4) * np.sin(y[i] * 0.4)
        veloc1[i,j] = veloc[i,j] * (1.0 + term * 0.1)
plt.contourf(veloc1)
plt.show()

f = open("veloctrue.in","w")
for i in range(nlat):
    for j in range(nlon):
        f.write("%f\n"%veloc1[i,j])
f.close()