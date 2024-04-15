#ifndef _SPMST_CSR_MATRIX_H
#define _SPMST_CSR_MATRIX_H

struct LSQRDict{
    
    float  atol,btol;
    float  conlim; 
    int itnlim,istop,itn;
    float  anorm,acond,arnorm;
    float  xnorm,rnorm;
    float  damp,weight;
    bool verbose;
    int nprocs;
    
    // some parameters defined
    // you could change it according to your own settings
    LSQRDict(int iterlim,float damp0,float smooth,int nproc_in = 1){
        atol = 1.0e-5,btol = 1.0e-5;
        conlim = 1.0e6;
        istop = 0;

        anorm = 0.0,acond = 0.0,arnorm = 0.0;
        xnorm = 0.0,rnorm = 0.0;

        itnlim = iterlim;
        damp = damp0;
        weight = smooth;

        verbose = false;

        nprocs = nproc_in;
    }
};

/**
 * Compressed Sparse Row matrix Class
 * Parameters:
 * ---------------------------------------
 *  for i-th row, the col number for nonzeros is indices[indptr[i]:indptr[i+1]]
 *  See Scipy doc for more information
 *  https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
 */
class csr_matrix{
private:
    int MATRIX_ROW,MATRIX_COL;
    bool ALLOCATE;

public:
    // row_index :indices
    int *indices,*indptr;
    float *data;
    int nonzeros;

    void initialize(int rows,int cols,int nar);

    // constructor
    csr_matrix(int rows,int cols,int nar);
    csr_matrix();

    // destructor
    ~csr_matrix();

    void lsqr_solver(const float* restrict b, float * restrict x,LSQRDict &dict) const;

    void read_binary(const char *filename);
    void write_binary(const char *filename) const;

    int rows() const;
    int cols() const;

    // merge function
    static void merge_csr_files(int nprocs,const char *outfile);
}; 

#endif