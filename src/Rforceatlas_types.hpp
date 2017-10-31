// Rforceatlas: Rcpp implementation of the ForceAtlas2 algorithm
// Copyright (C) 2017 José Tomás Atria jtatria at nomoi dot org
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
