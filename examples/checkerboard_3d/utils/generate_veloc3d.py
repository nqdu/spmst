import numpy as np
import matplotlib.pyplot as plt 

# input
xmin,xmax = 99.7,110.2
ymin,ymax = 25.7,35.3
nlon,nlat = 33,36
nlevelx = 0.6
nlevely = 0.6
nlevelz = 0.6
zstr = "0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6 6.4 6.8 7.2 7.6 8 8.4 8.8 9.2 9.6 15 35"
use_sph = False
shift_depth = False

# write velocity
z = list(map(lambda x: float(x),zstr.split()))
nz = len(z)
v = np.zeros((nz,nlat,nlon))
for i in range(nz):
    v[i,:,:] = 2.8 + 0.02 * z[i]

f = open("velocinit.in","w")
f.write("%d %d %d\n" %(nlon,nlat,nz))
f.write("%f %f\n%f %f\n" %(xmin,xmax,ymin,ymax))
f.write("%d %d\n" %(use_sph,shift_depth))
for i in range(nz):
    f.write("%f " %(z[i]))
f.write("\n")
for i in range(nz):
    for j in range(nlat):
        for k in range(nlon):
            f.write("%f\n"%(v[i,j,k]))
f.close()

xx = np.arange(nlon-2) + 1
yy = np.arange(nlat-2) + 1
zz = np.arange(nz-2) + 1

vtrue = v.copy()
for i in range(nz-2):
    for j in range(nlat-2):
        for k in range(nlon-2):
            v0 = 0.1 * np.sin(xx[k] * nlevelx) * np.sin(yy[j] * nlevely) * np.sin(zz[i] * nlevelz)
            vtrue[i+1,j+1,k+1] = v[i+1,j+1,k+1] * (1.0 + v0)

# plot 
plt.figure(1,figsize=(14,7))
plt.subplot(1,2,1)
plt.contourf(np.linspace(xmin,xmax,nlon),-np.array(z),(vtrue[:,3,:] - v[:,3,:]) / v[:,3,:] * 100)
plt.colorbar(orientation='horizontal')
plt.ylim(-10,0)
plt.subplot(1,2,2)
plt.contourf(np.linspace(xmin,xmax,nlon),np.linspace(ymin,ymax,nlat),(vtrue[3,:,:] - v[3,:,:]) / v[3,:,:] * 100)
plt.colorbar(orientation='horizontal')
plt.savefig("veloc3d.jpg")

f = open("veloctrue.in","w")
f.write("%d %d\n" %(nlon,nlat))
f.write("%f %f\n%f %f\n" %(xmin,xmax,ymin,ymax))
f.write("%d %d\n" %(use_sph,shift_depth))
for i in range(nz):
    f.write("%f " %(z[i]))
f.write("\n")
for i in range(nz):
    for j in range(nlat):
        for k in range(nlon):
            f.write("%f\n"%(vtrue[i,j,k]))
f.close()