#ifndef _WRK_H
#define _WRK_H 1

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "eigen_types.hpp"

using F = std::function<Vec(Vec,ind)>;

namespace wrk {
    void parallel( const Mat&, Mat&, const F );
    void serial( const Mat&, Mat&, const F );
};

#endif
