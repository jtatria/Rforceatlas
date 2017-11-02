#ifndef _EIGEN_TYPES
#define _EIGEN_TYPES 1

#include <RcppEigen.h>
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Mat;
typedef Eigen::Map<Eigen::SparseMatrix<double> > SpMat;
typedef SpMat::InnerIterator SpInIt;
typedef Eigen::VectorXd Vec;
typedef Mat::Scalar scalar;
typedef Eigen::Index ind;

#endif
