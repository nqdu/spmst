#include "shared/csr_matrix.hpp" 
#include <iostream>
#include <fstream>
#include "clsqr2/lsqr.hpp"
#include "shared/IO.hpp"
#include <string>

void csr_matrix:: 
initialize(int rows,int cols,int nar)
{
    MATRIX_ROW = rows;
    MATRIX_COL = cols;
    nonzeros = nar;
    indices = new int [nar]();
    data = new float [nar]();
    indptr = new int[rows + 1]();
    ALLOCATE = true;
}

csr_matrix::
csr_matrix(int rows,int cols,int nar){
    initialize(rows,cols,nar);
}
csr_matrix::csr_matrix(){ALLOCATE = false;}

// destructor
csr_matrix::
~csr_matrix()
{
    if(ALLOCATE){
        delete [] data;
        delete [] indptr;
        delete [] indices;
        ALLOCATE = false;
    }
}


void csr_matrix:: 
read_binary(const char *filename)
{
    // open file
    std::ifstream fp(filename,std::ios::binary);
    if (!fp.is_open()) {
        printf("Error opening file %s\n",filename);
        exit(1);
    }

    // scan first to read size of the matrix
    int m,n;
    fp.read((char*)&m,sizeof(int));
    fp.read((char*)&n,sizeof(int));

    // read nonzeros/rows in this block
    int idx,nar,nar1;
    nar1 = 0;
    while(!fp.read((char*)&idx,sizeof(int))){
        fp.read((char*)&nar,sizeof(int));
        nar1 += nar;
        
        // read temporay 
        float tmp[nar * 2];
        fp.read((char*)tmp,sizeof(float)*nar*2);
    }
    fp.close();

    // allocate space 
    this -> initialize(m,n,nar1);

    // now read matrix data
    fp.open(filename,std::ios::binary);
    fp.read((char*)&m,sizeof(int));
    fp.read((char*)&n,sizeof(int));
    while(!fp.read((char*)&idx,sizeof(int))){
        fp.read((char*)&nar,sizeof(int));
        int start = indptr[idx];
        indptr[idx + 1] = nar + start;
        fp.read((char*)(indices + start),sizeof(int)*nar);
        fp.read((char*)(data + start),sizeof(float)*nar);
    }

    // close file
    fp.close();
}

void csr_matrix:: 
write_binary(const char *filename) const
{
    // open file
    std::ofstream fp(filename,std::ios::binary);
    if(!fp.is_open()){
        printf("cannot open file %s\n",filename);
        exit(1);
    } 

    for(int i = 0; i < MATRIX_ROW; i++){
        int start = indptr[i], end = indptr[i+1];
        int nonzeros = end - start;
        fp.write((char*)&i,sizeof(int));
        fp.write((char*)&nonzeros,sizeof(int));
        fp.write((char*)(indices + start),sizeof(int)*nonzeros);
        fp.write((char*)(data + start),sizeof(int)*nonzeros);
    }

    fp.close();
}


/**
 * @brief solve sparse linear system Ax = b by lsqr method
 * 
 * @param b residual vector
 * @param x unknowns
 * @param dict LSQR dict
 */
void csr_matrix::
lsqr_solver(const float* b, float* __restrict x,LSQRDict &dict) const
{
    // initialize x 
    for(int i = 0; i < MATRIX_COL; i++) x[i] = 0.0;
    lsqr(MATRIX_ROW,MATRIX_COL,data,indices,indptr,dict.damp,b,x,dict.atol,
         dict.btol,dict.conlim,dict.itnlim,NULL,&dict.istop,&dict.itn,
         &dict.anorm,&dict.acond,&dict.rnorm,&dict.arnorm,&dict.xnorm,dict.nprocs);
}

int csr_matrix:: rows() const
{
    return MATRIX_ROW;
}

int csr_matrix:: cols() const
{
    return MATRIX_COL;
}

/**
 * merge csr files saved by several procs
 * @param nprocs # of threads used 
 * @param outfile # base filename, such as csr.bin
*/
void csr_matrix:: merge_csr_files(int nprocs,const char *outfile)
{
    std::ifstream fp;
    std::ofstream fpout;
    fpout.open(outfile,std::ios::binary);

    for(int irank = 0; irank < nprocs; irank ++){
        // open filename
        std::string filename = std::string(outfile) + "." + std::to_string(irank);
        fp.open(filename.c_str(),std::ios::binary);

        //get size 
        fp.seekg(0,std::ios::end);
        long size = (size_t)fp.tellg() - sizeof(int) * 2;
        fp.close();

        // read rows/cols 
        fp.open(filename.c_str(),std::ios::binary);
        int m,n;
        fp.read((char*)&m,sizeof(int));
        fp.read((char*)&n,sizeof(int));
        if(irank == 0) {
            fpout.write((char*)&m,sizeof(int));
            fpout.write((char*)&n,sizeof(int));
        }

        // write to fpout
        char *chunk = new char[size + 1];
        fp.read(chunk,size);
        fpout.write(chunk,size);
        fp.close();
        
        // remove this file
        std::remove(filename.c_str());
        delete [] chunk;
    }
    fpout.close();
}
