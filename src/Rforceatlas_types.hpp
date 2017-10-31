#ifndef TYPES_
#define TYPES_ 1

// [[Rcpp::plugins("cpp11")]]
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

typedef Rcpp::NumericMatrix RMatD;
typedef Rcpp::NumericVector RVecD;
typedef Rcpp::IntegerMatrix RMatI;
typedef Rcpp::IntegerVector RVecI;

typedef Eigen::Index ind;
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Mat::Scalar scalar;

#endif // TYPES_
