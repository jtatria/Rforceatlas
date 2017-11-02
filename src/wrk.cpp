// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "wrk.hpp"

namespace wrk {
    struct LambdaWorker : public RcppParallel::Worker {
        const Mat src;
        const F func;
        Mat tgt;

        LambdaWorker( const Mat& src, Mat& tgt, const F func )
            : src( src ), tgt( tgt ), func( func ) {}

        void operator()( std::size_t begin, std::size_t end ) {
            for( std::size_t i = begin; i < end; i++ ) {
                tgt.row( i ) = func( src.row( i ), i );
            }
        }
    };

    void parallel( const Mat& src, Mat& tgt, const F func ) {
        LambdaWorker wrkr( src, tgt, func );
        RcppParallel::parallelFor( 0, src.rows(), wrkr );
    }

    void serial( const Mat& src, Mat& tgt, const F func ) {
        for( ind i = 0; i < src.rows(); i++ ) {
            tgt.row( i ) = func( src.row( i ), i );
        }
    }
};
