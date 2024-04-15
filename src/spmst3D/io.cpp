#include "spmst3D/spmst3D.hpp"
#include "shared/IO.hpp"

/**
 * @brief write travel time to file
 * 
 * @param tsyn travel time data vector
 * @param filename output file
 */
void SPMST3D:: 
write_traveltime(const fvec &tsyn,const char *filename) const 
{
    FILE *fp = fopen(filename,"w");
    if(fp == NULL){
        printf("cannot open %s\n",filename);
        exit(1);
    }
    int nevents = stapairs.size();
    for(int ievt = 0; ievt < nevents; ievt ++){
        auto const &pair = stapairs[ievt];
        int nsta = pair.nreceivers;
        for(int ir = 0; ir < nsta; ir ++){
            float dist = compute_distance(pair.evlon,pair.evlat,
                                    pair.stlon[ir],pair.stlat[ir]);
            fprintf(fp,"%g %g %g %g %g %g %s\n",pair.evlon,pair.evlat,
                    pair.stlon[ir],pair.stlat[ir], dist,tsyn[pair.counter+ir],
                    pair.swdtp.c_str());
        }
    }
    fclose(fp);
}

/**
 * @brief read velocity model 
 * 
 * @param filename 
 * @param veloc_in shape (nz,nlat,nlon)
 * @param init_mesh if true, initialize spm mesh
 */
void SPMST3D:: 
read_model(const char *filename,fmat3 &veloc_in,bool init_mesh)
{
    std::ifstream fp; fp.open(filename);
    if(!fp.is_open()){
        printf("cannot open %s\n",filename);
        exit(1);
    }

    // read min/max lon/lat
    float lonmin,lonmax,latmin,latmax;
    read_file_param(fp,nlon,nlat,nz);
    read_file_param(fp,lonmin,lonmax);
    read_file_param(fp,latmin,latmax);
    read_file_param(fp,is_spherical,shift_depth);

    // create lon/lat vector
    model_lat.resize(nlat); model_lon.resize(nlon);
    for(int i = 0; i < nlon; i ++) {
        model_lon[i] = lonmin + (lonmax - lonmin) / (nlon-1.) * i;
    }
    for(int i = 0; i < nlat; i ++) {
        model_lat[i] = latmin + (latmax - latmin) / (nlat-1.) * i;
    }

    // read depth vector
    std::string line;
    std::getline(fp,line);
    std::vector<float> z(nz);
    std::vector<char> tmp; tmp.resize(line.size() + 1);
    memcpy(tmp.data(),line.data(),sizeof(char) * (line.size() + 1));
    char *starp = tmp.data(),*endp = NULL;
    for(int i = 0; i < nz; i ++) {
        z[i] = std::strtof(starp,&endp);
        starp = endp;
    }

    // print info 
    if(init_mesh) {
        printf("\nModel Description:\n");
        printf("===================================\n");
        if(is_spherical) {
            printf("using spherical coordinates ...\n");
            printf("lonmin = %f lonmax = %f\n",lonmin,lonmax);
            printf("latmin = %f latmax = %f\n",latmin,latmax);
            printf("nlat = %d nlon = %d nz = %d\n",nlat,nlon,nz);
        }
        else {
            printf("using cartesian coordinates ...\n");
            printf("xmin = %f xmax = %f\n",lonmin,lonmax);
            printf("ymin = %f ymax = %f\n",latmin,latmax);
            printf("ny = %d nx = %d nz = %d\n",nlat,nlon,nz);
        }
        printf("Grid points in depth direction:(km):\n");
        for(int i = 0; i < nz; i ++) {
            printf("%7.3f ",z[i]);
        }
        printf("\n");
        if(shift_depth) {
            printf("shifting depth and keep thickness ...\n");
        }
        else{
            printf("thickening/thinning the first layer ...\n");
        }

        // copy z to depth
        depth.resize(nz,nlat,nlon);
        for(int iz = 0; iz < nz; iz ++) {
        for(int iy = 0; iy < nlat; iy ++) {
        for(int ix = 0; ix < nlon; ix ++) {
            depth(iz,iy,ix) = z[iz];
        }}}
    }

    // resize depth and veloc_in
    veloc_in.resize(nz,nlat,nlon);
    for(int iz = 0; iz < nz; iz++){
    for(int ilat = 0; ilat < nlat; ilat ++){
    for(int ilon = 0; ilon < nlon; ilon ++){
        read_file_param(fp,veloc_in(iz,ilat,ilon));
    }}}

    // close file 
    fp.close();

    // initialize SPM Mesh
    if (init_mesh) {
        const int nrefine = 2;
        int nlonr = (nlon - 1) * nrefine;
        int nlatr = (nlat - 1) * nrefine;
        float dlonr = (lonmax - lonmin) / nlonr;
        float dlatr = (latmax - latmin) / nlatr;
        mesh.initialize(lonmin,latmin,dlonr,dlatr,nlonr,nlatr,is_spherical);
        mesh.create_graph(false);
    }
}


/**
 * @brief read observed data for SPMST  
 * 
 * @param filename observed data file
 */
void SPMST3D:: 
read_obsdata(const char *filename)
{
    printf("\nReading dispersion data:\n");
    printf("===================================\n");
    int ntRc,ntRg,ntLc,ntLg;
    char *starp,*endp;

    // open file
    std::string line;
    std::ifstream fp; fp.open(filename);
    if (!fp.is_open()) {
        printf("cannot open file %s\n",filename);
        exit(1);
    }

    // Rayleigh phase
    std::getline(fp,line);
    sscanf(line.c_str(),"%d",&ntRc);
    if(ntRc > 0) {
        printf("Rayleigh wave phase velocity used,periods:(s)\n");
        tRc.resize(ntRc);
        std::getline(fp,line); starp = (char*)line.data(), endp = NULL;
        for(int it = 0; it < ntRc; it ++) {
            tRc[it] = strtod(starp,&endp);
            starp = endp;
            printf("%g ",tRc[it]);
        }
        printf("\n");
    }

    // Rayleigh group
    std::getline(fp,line);
    sscanf(line.c_str(),"%d",&ntRg);
    if(ntRg > 0) {
        printf("Rayleigh wave group velocity used,periods:(s)\n");
        tRg.resize(ntRg);
        std::getline(fp,line); starp = (char*)line.data(), endp = NULL;
        for(int it = 0; it < ntRg; it ++) {
            tRg[it] = strtod(starp,&endp);
            starp = endp;
            printf("%g ",tRg[it]);
        }
        printf("\n");
    }

    // Love phase
    std::getline(fp,line);
    sscanf(line.c_str(),"%d",&ntLc);
    if(ntLc > 0) {
        printf("Love wave phase velocity used,periods:(s)\n");
        std::getline(fp,line); starp = (char*)line.data(), endp = NULL;
        tLc.resize(ntLc);
        for(int it = 0; it < ntLc; it ++) {
            tLc[it] = strtod(starp,&endp);
            starp = endp;
            printf("%g ",tLc[it]);
        }
        printf("\n");
    }

    // Love group
    std::getline(fp,line);
    sscanf(line.c_str(),"%d",&ntLg);
    if(ntLg > 0) {
        printf("Love wave group velocity used,periods:(s)\n");
        tLg.resize(ntLg);
        std::getline(fp,line); starp = (char*)line.data(), endp = NULL;
        for(int it = 0; it < ntLg; it ++) {
            tLg[it] = strtod(starp,&endp);
            starp = endp;
            printf("%g ",tLg[it]);
        }
        printf("\n");
    }

    // some dummy variables
    char dummy;
    float x,y; 
    int nsta,pid,waveid,veltp;

    // scan the file for the first time to find how many station pairs 
    int nevts = 0;
    while(std::getline(fp,line)){
        sscanf(line.c_str(),"%c%f%f%d",&dummy,&x,&y,&nsta);

        // skip all receivers
        for(int i = 0; i < nsta; i++){
            std::getline(fp,line);
        }
        nevts += 1;
    }
    fp.close();
    stapairs.reserve(nevts);

    // scan for the second time to read all the information
    int ndata = 0;
    fp.open(filename);
    std::getline(fp,line); if (ntRc > 0) std::getline(fp,line); 
    std::getline(fp,line); if (ntRg > 0) std::getline(fp,line); 
    std::getline(fp,line); if (ntLc > 0) std::getline(fp,line); 
    std::getline(fp,line); if (ntLg > 0) std::getline(fp,line); 
    while(std::getline(fp,line)){
        StationPair pair;
        sscanf(line.c_str(),"%c%f%f%d%d%d%d",&dummy,&x,&y,&nsta,&pid,&waveid,&veltp);
        pair.evlat = y; pair.evlon = x;
        pair.period_id = pid;
        pair.nreceivers = nsta;
        pair.stlon.resize(nsta); pair.stlat.resize(nsta); 
        pair.ttime.resize(nsta);
        pair.counter = ndata;

        // get wave type
        std::string wtp;
        if ( waveid == 2 && veltp ==0 ) 
            wtp = "Rc";
        else if ( waveid == 2 && veltp==1 ) 
            wtp = "Rg";
        else if ( waveid == 1 && veltp==0 ) 
            wtp = "Lc";
        else {
            wtp = "Lg";
        }
        pair.swdtp = wtp;

        // read all receivers
        for(int i = 0; i < nsta; i++){
            std::getline(fp,line);
            float vel;
            sscanf(line.c_str(),"%f%f%f",&pair.stlon[i],&pair.stlat[i],&vel);
            float dist = this -> compute_distance(pair.evlon,pair.evlat,
                                                pair.stlon[i],pair.stlat[i]);
            pair.ttime[i] = dist / vel; 
        }

        // append to stapairs
        stapairs.push_back(pair);

        // renew ndata
        ndata += nsta;
    }
    fp.close();

    // allocate time vector space
    tobs.resize(ndata);
    for(int ievt = 0; ievt < nevts; ievt ++){
        int c = stapairs[ievt].counter;
        int nsta = stapairs[ievt].nreceivers;
        memcpy(&tobs[c],stapairs[ievt].ttime.data(),sizeof(float)*nsta);
    }
    printf("# of observations = %d\n",ndata);
}

/**
 * @brief save current model to txt file
 * 
 * @param filename 
 * @param veloc shape(nz,nlat,nlon)
 */
void SPMST3D:: 
save_model(const char *filename,const fmat3 &veloc) const
{
    FILE *fp = fopen(filename,"w");
    for(int iz = 0; iz < nz; iz ++) {
    for(int iy = 0; iy < nlat; iy ++) {
    for(int ix = 0; ix < nlon; ix ++) {
        fprintf(fp,"%g %g %g %g\n",model_lon[ix],model_lat[iy],depth(iz,iy,ix),veloc(iz,iy,ix));
    }}}
    fclose(fp);
}