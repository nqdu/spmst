import numpy as np
import matplotlib.pyplot as plt 

# input
xmin,xmax = 89.2,  90.3
ymin,ymax = 38.9,   39.3
nlon,nlat = 40,40
nlevelx = 0.5
nlevely = 0.5
use_sph = True

# write velocity
v = np.zeros((nlat,nlon))
v += 3.0

f = open("velocinit.in","w")
f.write("%d %d\n" %(nlon,nlat))
f.write("%f %f\n%f %f\n" %(xmin,xmax,ymin,ymax))
f.write("%d\n" %(use_sph))
for j in range(nlat):
    for k in range(nlon):
        f.write("%f\n"%(v[j,k]))
f.close()

xx = np.arange(nlon-2) + 1
yy = np.arange(nlat-2) + 1

vtrue = v.copy()
for j in range(nlat-2):
    for k in range(nlon-2):
        v0 = 0.1 * np.sin(xx[k] * nlevelx) * np.sin(yy[j] * nlevely)
        vtrue[j+1,k+1] = v[j+1,k+1] * (1.0 + v0)
f.close()

# plot 
plt.contourf(np.linspace(xmin,xmax,nlon),np.linspace(ymin,ymax,nlat),(vtrue[:,:] - v[:,:]) / v[:,:] * 100)
plt.colorbar(orientation='horizontal')
plt.savefig("veloc2d.jpg")

f = open("veloctrue.in","w")
f.write("%d %d\n" %(nlon,nlat))
f.write("%f %f\n%f %f\n" %(xmin,xmax,ymin,ymax))
f.write("%d\n" %(use_sph))
for j in range(nlat):
    for k in range(nlon):
        f.write("%f\n"%(vtrue[j,k]))
f.close()