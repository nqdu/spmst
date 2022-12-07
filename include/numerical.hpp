#ifndef _SPM2D_NUMERICAL_H
#define _SPM2D_NUMERICAL_H

#include <limits>
#include <unsupported/Eigen/CXX11/Tensor>

const float inf = std::numeric_limits<float>::infinity();
typedef Eigen::Array<float,-1,-1,Eigen::RowMajor> fmat2;
typedef Eigen::Array<int,-1,-1,Eigen::RowMajor> imat2;
typedef Eigen::Array<float,-1,1> fvec;
typedef Eigen::Array<int,-1,1> ivec;

#endif