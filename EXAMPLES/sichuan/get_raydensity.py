import numpy as np
import sys 

def get_great_circle(lon0,lat0,lon1,lat1,nseg):
    import pyproj 
    g = pyproj.Geod(a=1,f=0)
    data = np.zeros((nseg,2))
    data[0,:] = [lon0,lat0]
    data[1:nseg-1,:] = np.array(g.npts(lon0,lat0,lon1,lat1,nseg-2))
    data[-1,:] = [lon1,lat1]
    return data 


def main():
    if len(sys.argv) !=6:
        print("usage: ./this xmin xmax ymin ymax is_spherical")
        print(sys.argv)
        exit(1)
    xmin,xmax,ymin,ymax = map(lambda x: float(x),sys.argv[1:5])
    is_spherical = int(sys.argv[5])
    nx = 128
    ny = 128
    x = np.linspace(xmin,xmax,nx)
    y = np.linspace(ymin,ymax,ny)
    dx = x[1] - x[0]; dy = y[1] - y[0]
    ray = np.zeros((ny,nx)) + 0.
    nseg = 100

    # now loop the data file to get all density 
    f = open("surfdata.txt","r")
    line = f.readline()
    while line:
        info = line.split()
        evlo,evla = map(lambda x:float(x),info[1:3])
        nsta = int(info[3])
        for _ in range(nsta):
            line = f.readline()
            info = line.split()
            stlo,stla = map(lambda x:float(x),info[:2])
            data = get_great_circle(evlo,evla,stlo,stla,nseg)
            for i in range(nseg):
                x0,y0 = data[i,:]
                ix = int((x0-xmin)/dx)
                iy = int((y0-ymin) / dy)
                ray[iy:iy+2,ix:ix+2] += 1.0
        
        line = f.readline()
    f.close()

    ray[ray > 0.0] = 1.0
    ray[ray == 0.0] = np.nan 

    # write out ray
    f = open("white.dat","w")
    for i in range(ny):
        for j in range(nx):
            f.write("%f %f %f\n"%(x[j],y[i],ray[i,j]))
    f.close()

main()

        