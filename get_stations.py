import numpy as np 
from obspy.geodetics.base import locations2degrees

nsta = 45 
x = 89.21 + np.random.rand(nsta) * (90.29 - 89.21)
y = 38.91 + np.random.rand(nsta) * (39.29 - 38.91)

f = open("station.txt","w")

for i in range(nsta-1):
    f.write("# %f %f %d\n" %(x[i],y[i],nsta-1-i))
    for j in range(i+1,nsta):
        dist = locations2degrees(y[i],x[i],y[j],x[j]) * 6371.0 * np.pi / 180.
        f.write("%f %f %f\n" %(x[j],y[j],dist / 3.0))
f.close()


