import sys 
import numpy as np
from scipy.interpolate import griddata
import os 
os.environ["MKL_NUM_THREADS"] = "1" 

def lonlat2xyz(lon,lat,dep):
    r = 6371.
    r2d = np.pi / 180
    x = (r - dep) * np.cos(lat*r2d) * np.cos(lon*r2d)
    y = (r - dep) * np.cos(lat*r2d) * np.sin(lon*r2d)
    z = (r - dep) * np.sin(lat*r2d)

    return x,y,z

def get_xyz(data,is_spherical):
    if not is_spherical:
        x = data[:,0]
        y = data[:,1]
        z = data[:,2]

        return x,y,z 
    else:
        x,y,z = lonlat2xyz(data[:,0],data[:,1],data[:,2])

        return x,y,z

def load_topo(topo_file):
    topo = np.loadtxt(topo_file,skiprows=2)
    f = open(topo_file,"r")
    line = f.readline()
    nz,nx = map(lambda x:int(x),line.split())
    line = f.readline()
    xmin,xmax,zmin,zmax = map(lambda x:float(x),line.split())
    n = nz * nx

    xtopo = np.zeros((n))
    ztopo = np.zeros((n))

    for i in range(nz):
        for j in range(nx):
            id = i * nx + j 
            xtopo[id] = xmin + (xmax - xmin) / (nx-1) * j
            ztopo[id] = zmin + (zmax - zmin) / (nz-1) * i
    
    topo *= -0.001
    return xtopo,ztopo,topo

def main():
    if len(sys.argv) != 4:
        print("Usage: ./this z0 model_name topo_file")
        exit(1)
    z0 = float(sys.argv[1])
    model_name = sys.argv[2]
    topo_file = sys.argv[3]
    
    # get spherical flag
    is_spherical = 0
    f = open("spmtomo.in","r")
    lines = f.readlines()
    f.close()
    n0 = 0
    for line in lines:
        if line[0] == '#': continue
        if n0 !=3:
            n0 += 1
        else:
            is_spherical = int(line.split()[0])
            break 

    # load model
    model = np.loadtxt(model_name)
    model_x,model_y,model_z = get_xyz(model[:,:3],is_spherical)


    # load topography
    topo_x,topo_y,topo_z = load_topo(topo_file)

    # profile
    n = 128

    # get coordinate
    cord = np.zeros((n*n,3))
    lon0 = model[:,0].min(); lon1 = model[:,0].max()
    lat0 = model[:,1].min(); lat1 = model[:,1].max()
    for i in range(n):
        lat = lat0 + (lat1 - lat0) / (n-1) * i
        for j in range(n):
            lon = lon0 + (lon1 - lon0) / (n-1) * j
            id = i * n + j
            cord[id,0] = lon
            cord[id,1] = lat
            cord[id,2] = z0

    # get topography at cord
    loc_topo = griddata((topo_x,topo_y),topo_z,(cord[:,0],cord[:,1]))
    idx = np.where(np.isnan(loc_topo))[0]
    loc_topo[idx] = griddata((topo_x,topo_y),topo_z,(cord[idx,0],cord[idx,1]),'nearest')

    cord[:,2] += loc_topo 
    print(np.mean(cord[:,2]))
 
    # get velocity
    cordx,cordy,cordz = get_xyz(cord,is_spherical)
    v = griddata((model_x,model_y,model_z),model[:,-1],(cordx,cordy,cordz))
    idx = np.where(np.isnan(v))[0]
    v[idx] = griddata((model_x,model_y,model_z),model[:,-1],(cordx[idx],cordy[idx],cordz[idx]))

    f = open("out_v.txt","w")
    for i in range(n):
        for j in range(n):
            id = i * n + j
            f.write("%f %f %f\n"%(cord[id,0],cord[id,1],v[id]))
    f.close()

main()