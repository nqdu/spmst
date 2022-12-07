#include "shared/csr_matrix.hpp" 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "clsqr2/lsqr.hpp"

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

// reader
void csr_matrix::
read(const char *filename)
{
    // open file
    FILE *fp;
    if((fp=fopen(filename,"r"))==NULL){
        printf("cannot open file %s\n",filename);
        exit(0);
    } 

    // read from file
    char line[100],dummy;
    int rw_idx,nar;
    
    while(fgets(line,sizeof(line),fp)!=NULL){
        if(line[0] != '#') break;
        sscanf(line,"%c%d%d",&dummy,&rw_idx,&nar);
        int start = indptr[rw_idx];
        indptr[rw_idx + 1] = nar + start;
        int end = indptr[rw_idx + 1];
        //std::cout << start << " " << end << std::endl;
        for(int i=start;i<end;i++){
            assert(fscanf(fp,"%d%f\n",indices +i ,data + i) == 2);
        }
    }
}

void csr_matrix::
lsqr_solver(const float* restrict b, float* restrict x,LSQRDict &dict)
{
    // initialize x 
    for(int i = 0; i < MATRIX_COL; i++) x[i] = 0.0;
    lsqr(MATRIX_ROW,MATRIX_COL,data,indices,indptr,dict.damp,b,x,dict.atol,
         dict.btol,dict.conlim,dict.itnlim,NULL,&dict.istop,&dict.itn,
         &dict.anorm,&dict.acond,&dict.rnorm,&dict.arnorm,&dict.xnorm);
}

int csr_matrix:: rows(){
    return MATRIX_ROW;
}

int csr_matrix:: cols(){
    return MATRIX_COL;
}