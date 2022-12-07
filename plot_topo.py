import numpy as np
import matplotlib.pyplot as plt 

# load topo
nlon = 265; nlat = 97
topo = np.loadtxt("topo.dat",skiprows=2).reshape((nlat,nlon)) * 0.001
lon = np.linspace(89.2,90.3,nlon); lat = np.linspace(38.9,39.3,nlat)
X,Y = np.meshgrid(lon,lat)
print(X.shape,Y.shape)
# plot raypath
ray = np.loadtxt("test.dat")

# plot
ax = plt.figure().add_subplot(projection='3d')
ax.plot_surface(X,Y,topo)

ax.scatter(ray[:,0],ray[:,1],ray[:,2],color='r')
plt.show()