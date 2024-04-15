import numpy as np
import sys 
import os 
import matplotlib.pyplot as plt 

def get_topo(H,R):
    nlat = 51 
    nlon = 51 
    lat = np.linspace(-5,5,nlat)
    lon = np.linspace(-5,5,nlon)

    topo = np.zeros((nlat,nlon))
    for i in range(nlat):
        for j in range(nlon):
            t = (lat[i] )**2 + (lon[j])**2
            topo[i,j] = H * np.exp(- t / R**2 )


    f = open("topo.txt","w")
    f.write(f"{nlat} {nlon}\n")
    f.write("-5 5 -5 5\n")
    for i in range(nlat):
        for j in range(nlon):
            f.write("%f\n"%(topo[i,j]))
    f.close()

def main():

    R = np.linspace(0.5,1.5,11)
    H = np.linspace(200,1000,17)
    nh = len(H)
    nr = len(R)
    ratio = np.zeros((nr,nh))

    for i in range(nr):
        print(i)
        for j in range(nh):
            get_topo(H[j],R[i])
            os.system("../../bin/syn spmst.in surfdata.txt topo.txt > out.log")

            d = np.loadtxt("ray.dat",comments='>')
            dist0 = np.sqrt(np.sum((d[-1,:] - d[0,:])**2))
            dist1 = np.sum(np.sqrt(np.sum(np.diff(d,axis=0)**2,axis=1)))

            ratio[i,j] = (dist1 - dist0) / dist0 

    # write out
    f = open('ratio.dat',"w")
    for i in range(nr):
        for j in range(nh):
            f.write("%f %f %f\n"%(H[j],ratio[i,j],R[i]))
        f.write(">\n")
    f.close()

    for i in range(nr):
        plt.plot(H,ratio[i,:],label='%g'%(R[i]))
    
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()