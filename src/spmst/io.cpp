#include "spmst.hpp"
#include "shared/gps2dist.hpp"

void SPMST:: 
write_traveltime(const fvec &tsyn,const char *filename) const 
{
    FILE *fp = fopen(filename,"w");
    if(fp == NULL){
        printf("cannot open %s\n",filename);
        exit(1);
    }
    int inum = 0;
    for(int ievt = 0; ievt < nevents; ievt ++){
        int nsta = nrecvs_per_event[ievt];
        for(int ir = 0; ir < nsta; ir ++){
            float dist = gps2dist(evlon[ievt],stalon(ievt,ir),evlat[ievt],
                                 stalat(ievt,ir),earth);
            fprintf(fp,"%g %g %g %g %g %g\n",evlon[ievt],evlat[ievt],
                    stalon(ievt,ir),stalat(ievt,ir),
                    dist,dist/tsyn[inum+ir]);
        }
        inum += nsta;
    }
    fclose(fp);
}