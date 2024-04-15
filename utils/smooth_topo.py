from scipy.ndimage import gaussian_filter
import numpy as np
import matplotlib.pyplot as plt 
import sys

def load_topo(topo_file):
    topo = np.loadtxt(topo_file,skiprows=2)
    f = open(topo_file,"r")
    line = f.readline()
    nx,nz = map(lambda x:int(x),line.split()[:2])
    line = f.readline()
    xmin,xmax = map(lambda x:float(x),line.split()[:2])
    line = f.readline()
    zmin,zmax = map(lambda x:float(x),line.split()[:2])

    xtopo = np.zeros((nz,nx))
    ztopo = np.zeros((nz,nx))
    topo = np.reshape(topo,(nz,nx))

    for i in range(nz):
        for j in range(nx):
            xtopo[i,j] = xmin + (xmax - xmin) / (nx-1) * j
            ztopo[i,j] = zmin + (zmax - zmin) / (nz-1) * i
    
    topo *= -0.001
    return xtopo,ztopo,topo

def save_topo(xtopo,ztopo,topo,filename):
    xmin = xtopo.min()
    xmax = xtopo.max()
    zmin = ztopo.min()
    zmax = ztopo.max()
    nz,nx = xtopo.shape 

    f = open(filename,"w")
    f.write("%d %d\n"%(nx,nz))
    f.write('%f %f\n%f %f\n'%(xmin,xmax,zmin,zmax))

    for i in range(nz):
        for j in range(nx):
            f.write("%f\n"%(-topo[i,j]*1000))
    f.close()


def main():
    if len(sys.argv) != 2:
        print("Usage: python smooth_topy.py topofile")
        exit(1)
    
    # read topofile
    topofile = sys.argv[1]
    xtopo,ztopo,topo = load_topo(topofile)

    topo1 = gaussian_filter(topo,[2.5,1.5])

    plt.figure(1,figsize=(22,6))
    plt.subplot(1,2,1)
    plt.contourf(xtopo,ztopo,topo)
    plt.colorbar()

    plt.subplot(1,2,2)
    plt.contourf(xtopo,ztopo,topo1)
    plt.colorbar()
    plt.show()

    # save topo
    save_topo(xtopo,ztopo,topo1,topofile + ".smooth")