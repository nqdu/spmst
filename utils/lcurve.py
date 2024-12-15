import numpy as np
from numba import jit 
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import sys

@jit(nopython=True)
def add_smooth3D(nlon,nlat,nz,indices,indptr,data):
    # get parameters required
    n = nlon * nlat # model dimension
    nonzeros = len(data)
    m = 0

    count = 0
    nar = 0
    for k in range(nz):
        for i in range(nlat):
            for j in range(nlon):
                idx = k * nlon * nlat + i * nlon + j
                if j ==0 or j == nlon-1 or i==0 or i == nlat-1 or k ==0 or k == nz-1:
                    rwc = count + m
                    data[nar] = 4.0
                    indices[nar] = idx
                    indptr[rwc + 1] = indptr[rwc] + 1
                    nar += 1
                    count += 1
                    continue
                else:
                    rwc = count + m;  # current row
                    indptr[rwc +1] = indptr[rwc] + 7
                    clc = idx
                    data[nar] = 6.0
                    indices[nar] = clc

                    # x direction
                    data[nar + 1] = -1.0
                    indices[nar + 1] = clc - 1
                    data[nar + 2] = -1.0
                    indices[nar + 2] = clc + 1

                    # y direction
                    data[nar + 3] = -1.
                    indices[nar + 3] = clc - nlon
                    data[nar + 4] = -1.
                    indices[nar + 4] = clc + nlon

                    # z direction
                    data[nar + 5] = -1.
                    indices[nar + 5] = clc - nlon * nlat
                    data[nar + 6] = -1.
                    indices[nar + 6] = clc + nlon * nlat

                    nar += 7
                    count += 1
    nonzeros = nar

    return nonzeros


@jit(nopython=True)
def add_smooth2D(nlon,nlat,indices,indptr,data):
    # get parameters required
    n = nlon * nlat # model dimension
    nonzeros = len(data)
    m = 0

    count = 0
    nar = 0
    for j in range(nlat):
        for i in range(nlon):
            idx = j * nlon + i
            if j ==0 or j == nlat-1 or i==0 or i == nlon-1:
                rwc = count + m
                data[nar] = 4.0
                indices[nar] = idx
                indptr[rwc + 1] = indptr[rwc] + 1
                nar += 1
                count += 1
                continue
            else:
                rwc = count + m;  # current row
                indptr[rwc +1] = indptr[rwc] + 5
                clc = idx
                data[nar] = 4.0
                indices[nar] = clc

                # x direction
                data[nar + 1] = -1.0
                indices[nar + 1] = clc - 1
                data[nar + 2] = -1.0
                indices[nar + 2] = clc + 1

                # y direction
                data[nar + 3] = -1.
                indices[nar + 3] = clc - nlon
                data[nar + 4] = -1.
                indices[nar + 4] = clc + nlon

                nar += 5
                count += 1
    nonzeros = nar

    return nonzeros

@jit(nopython=True)
def convert(md,nx,ny,nz):
    n = nx * ny * (nz - 1)
    x = np.zeros((n),dtype=float)
    for i in range(nz - 1):
        for j in range(ny):
            for k in range(nx):
                idx = i * ny * nx + j * nx + k
                x[idx] = md[i*ny*nx+j*nx+k]
    return x

def get_laplacian3D(nx,ny,nz):
    # generate csr type laplacian matrix 

    n = nx * ny * (nz-1)
    nonzeros = n * 7
    indptr = np.zeros((n+1),dtype=int)
    indices = np.zeros((nonzeros),dtype=int )
    data = np.zeros((nonzeros),dtype = float)
    nonzeros = add_smooth3D(nx,ny,nz-1,indices,indptr,data)
    indices = indices[:nonzeros]
    data = data[:nonzeros]
    smat = csr_matrix((data, indices, indptr),shape=(n,n))

    return smat

def get_laplacian2D(nx,ny):
    # generate csr type laplacian matrix 

    n = nx * ny
    nonzeros = n * 5
    indptr = np.zeros((n+1),dtype=int)
    indices = np.zeros((nonzeros),dtype=int )
    data = np.zeros((nonzeros),dtype = float)
    nonzeros = add_smooth2D(nx,ny,indices,indptr,data)
    indices = indices[:nonzeros]
    data = data[:nonzeros]
    smat = csr_matrix((data, indices, indptr),shape=(n,n))

    return smat


def compute_roughness(nx,ny,nz,smat,result_dir:str):
    # load model and compute m - m0
    md = np.loadtxt(result_dir + '/mod_iter1.dat') [:,-1]
    md -= np.loadtxt(result_dir + '/mod_iter0.dat') [:,-1]
    x = convert(md,nx,ny,nz)
    roughness_smooth = np.sqrt(np.sum((smat * x)**2))
    roughness_damp = np.sqrt(np.sum((x)**2))

    return roughness_smooth,roughness_damp

def main():
    if len(sys.argv) == 1:
        print("compute |m-m0|, |L(m-m0)| for mod_iter1.dat")
        print("Usage: python lcurve.py init_model model_dir")
        print("example: python lcurve.py MOD.init model")

        exit(1)
    
    # get input args
    init_model_str = sys.argv[1]
    result_dir = sys.argv[2]

    # read model dimension/
    # read corner 
    f = open(init_model_str,"r")
    line = f.readline()
    info = line.split()
    if len(info) < 3:
        print("2-D case ...")
        nx,ny = map(lambda x: int(x), line.split()[:2])
        nz = 2
        smat = get_laplacian2D(nx,ny)
    else:
        print("3-D case ...")
        nx,ny,nz = map(lambda x: int(x), line.split()[:3])
        smat = get_laplacian3D(nx,ny,nz)
    f.close()


    # get roughness
    r_smooth,r_damp = compute_roughness(nx,ny,nz,smat,result_dir)

    print("roughness smooth |L(m-m0)| = %g" %(r_smooth))
    print("roughness damp |(m-m0)| = %g" %(r_damp))
    
main()