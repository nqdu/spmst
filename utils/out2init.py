import numpy as np
import sys 

def main():
    if len(sys.argv) != 3:
        print("Usage: out2init.py models/mod_iter10.dat velocinit.dat")
        exit(1)

    outfile = sys.argv[1]
    reffile = sys.argv[2]

    # read outfile
    d = np.loadtxt(outfile)[:,-1]

    f = open(reffile,"r")
    lines = f.readlines()
    f.close()

    f = open(reffile + '.out',"w")
    f.writelines(lines[:5])
    for i in range(d.shape[0]):
        f.write("%f\n"%d[i])
    f.close()

main()