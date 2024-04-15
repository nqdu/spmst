import numpy as np
import matplotlib.pyplot as plt 

def gaus2d_mountain(x=0, y=0, mx=0, my=0, sx=0.1, sy=0.1):
    return 1 * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))

nlon = 64
nlat = 64
v = np.zeros((nlat,nlon)) + 1.5
x = np.linspace(-5,5,nlon)
y = np.linspace(-5,5,nlat)
am = 0.8

for i in range(nlat):
    for j in range(nlon):
        dv = 0.2 * np.sin(am * j) * np.sin(am * i)
        dv1 = 20 * gaus2d_mountain(x[j]+3,y[i]+3)
        dv2 = 20 * gaus2d_mountain(x[j]+3,y[i]-3)
        #dv1 = 3.0 * (x[j]  0.123)
        #dv2 = 0.
        v[i,j] += dv1 + dv2

f = open("veloc.in","w")
f.write("%d %d\n"%(nlon,nlat))
for i in range(nlat):
    for j in range(nlon):
        f.write("%f %f %f\n"%(x[j],y[i],v[i,j]))
f.close()

plt.contourf(x,y,v)
plt.colorbar()
plt.show()
