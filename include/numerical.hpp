#ifndef _SPMST_NUMERICAL_H
#define _SPMST_NUMERICAL_H

#include <limits>

#if __GNUC__ > 11
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
    #include <unsupported/Eigen/CXX11/Tensor>
#pragma GCC diagnostic pop
#else 
    #include <unsupported/Eigen/CXX11/Tensor>
#endif


const float inf = std::numeric_limits<float>::infinity();
typedef Eigen::Array<float,-1,-1,Eigen::RowMajor> fmat2;
typedef Eigen::Tensor<float,3,Eigen::RowMajor> fmat3;
typedef Eigen::Tensor<double,3,Eigen::RowMajor> dmat3;
typedef Eigen::Array<int,-1,-1,Eigen::RowMajor> imat2;
typedef Eigen::Array<float,-1,1> fvec;
typedef Eigen::Array<double,-1,1> dvec;
typedef Eigen::Array<int,-1,1> ivec;

#endif