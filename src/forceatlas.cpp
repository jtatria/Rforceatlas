//' @backref src/nufa.cpp

#include "Rforceatlas_types.hpp"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace Rcpp;

using F = std::function<Vec(Vec,ind)>;

const scalar EPSILON = 0.00000001;

volatile static int WTF;
void wtf() {
    WTF = -99;
}

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

// [[Rcpp::export]]
scalar attr_func(
    const scalar k, const scalar d, const scalar deg_i, const scalar wgt,
    const bool linlog, const bool nohubs, const bool overlap
) {
    if( overlap ) return 0;
    if( wgt == 0 ) return 0;
    scalar a = k * wgt;
    a *= linlog ? std::log<scalar>( 1.0 + d ).real() / d : 1.0;
    a /= nohubs ? deg_i : 1.0;
    return a;
}

// [[Rcpp::export]]
scalar repl_func(
    const scalar k, const scalar d, const scalar deg_i, const scalar deg_j,
    const bool overlap
) {
    scalar r = ( overlap ) ? 100 * k * deg_i * deg_j : k * deg_i * deg_j / d / d;
    return r * -1.0;
}

// [[Rcpp::export]]
scalar grav_func(
    const scalar k, const scalar d, const scalar deg_i, const scalar G, const bool strong
) {
    scalar g = k * G * ( deg_i + 1.0 ) / d;//( strong ) ? k * deg_i * G : k * deg_i * G / d;
    return g;
}

bool get_overlap( const scalar d, const ind i, const ind j ) {
    return false;
}

scalar edge_weight( scalar delta, scalar w ) {
    if( w == 0.0 ) return 0.0;
    if( delta <= 0.0 ) return 1.0;
    if( delta <= 1.0 ) return w;
    return std::pow<scalar>( w, delta );
}

Mat force_func(
    const Mat& pos, const Vec& deg, const Mat& W, const scalar k,
    const scalar delta, const scalar G,
    const bool linlog, const bool nohubs, const bool strong,
    const Vec& orig
) {
    Mat force = Mat::Zero( pos.rows(), pos.cols() );
    // F inner_l = [&]( const Vec& vi, ind i ) -> Vec {
    for( ind i = 0; i < pos.rows(); i++ ) {
        Vec vi = pos.row( i );
        scalar deg_i = deg[i];
        for( ind j = 0; j < pos.rows(); j++ ) {
            if( i == j ) continue;
            Vec vj = pos.row( j );
            scalar deg_j = deg[j];
            Vec ij = ( vi - vj ).eval();
            scalar d = ij.norm();
            bool overlap = get_overlap( d, i, j );
            scalar w = edge_weight( W( i, j ), delta );
            scalar a = 0; //attr_func( k, d, deg_i, w, linlog, nohubs, overlap );
            scalar r = 0; //repl_func( k, d, deg_i, deg_j, overlap );
            Vec dj = ( ij * ( a + r ) );
            force.row( i ) += dj;
        }
        Vec vg = ( orig - vi );
        scalar d = vg.norm();
        scalar g = grav_func( k, d, deg_i, G / k, strong );
        Vec dg  = ( vg * g );
        force.row( i ) += dg;
    }
    //     return ( fi + gravity( vi, orig, G, degi, strong ) );
    // };
    // impl::serial( pos, force, inner_l );
    if( ( force.array() == 0 ).any() ) wtf();
    return force;
}

bool conv_theta( const Vec& swing ) {
    return false;
}

Vec speed_func(
    const Mat& f_t0, const Mat& f_t1, const Vec& deg,
    const scalar tol, const scalar ks, const scalar ksmax
) {
    // Vec swing = ( f_t0 - f_t1 ).rowwise().norm();
    // if( ( conv_theta( swing ) ) ) return Vec::Zero( f_t0.rows() );
    // Vec tract = ( f_t0 + f_t1 ).rowwise().norm() / 2;
    // scalar swing_g = swing.cwiseProduct( ( deg.array() + 1.0 ).matrix() ).sum();
    // scalar tract_g = tract.cwiseProduct( ( deg.array() + 1.0 ).matrix() ).sum();
    // scalar speed_g = tol * tract_g / swing_g;
    // Vec speed = ( swing.cwiseSqrt() * ( speed_g + 1.0 ) ).array().pow( -1 ) * ( ks * speed_g );
    // Vec speed_max = ( f_t0.rowwise().norm() ).array().pow( -1.0 ) * ksmax;
    // return ( speed.array() > speed_max.array() ).select( speed_max, speed );
    return Vec::Ones( f_t0.rows() );
}

void disp_func( const Mat& force, const Vec& speed, Mat& disp ) {
    for( ind i = 0; i < disp.rows(); i++ ) {
        disp.row( i ) = force.row( i ) * speed[i];
    }
    disp.eval();
}

//' ForceAtlas2 graph layout algortihm
//'
//' Place vertices on an n-dimensional space using the ForceAtlas2 algorithm by
//' Jacomy et al. (2014), originally developed in Java for the Gephi graph analysis software.
//' This implementation has been written from scratch in C++ following the R implementation by
//' Klockiewicz and Álvarez.
//'
//' A full description of the algorithm and all of its parameters can be found in the reference
//' below.
//'
//' @references
//' Jacomy, Mathieu, Tommaso Venturini, Sebastien Heymann, and Mathieu Bastian. "ForceAtlas2, a
//' Continuous Graph Layout Algorithm for Handy Network Visualization Designed for the Gephi
//' Software." PLOS ONE 9, no. 6 (June 10, 2014): e98679.
//' \link{https://doi.org/10.1371/journal.pone.0098679}.
//'
//' @param m      An adjacency matrix.
//' @param iter   Integer, number of iterations.
//' @param linlog Logical. Use log( distance ) for initial attraction forces instead of plain
//'               distance. Defaults to FALSE.
//' @param nohubs Logical. Dissuade hubs, placing authorities closer to the center. Defaults to
//'               FALSE.
//' @param k      Numeric. Repulsion scaling factor. Defaults to 400.
//' @param G      Numeric. Gravity scaling factor. Defaults to 1.
//' @param ks     Numeric. Local speed scaling factor. Defaults to .1.
//' @param ksmax  Numeric. Maximum local speed. Defaults to 10.
//' @param delta  Numeric. Exponent to weight attraction forces according to the entries in m.
//' @param tol    Numeric. Tolerance to swinging. Defaults to .1
//' @param dim    Integer. Number of dimensions to position vertices against. Defaults to 2.
//' @param init   Numeric matrix of dimensions nrow( m ) * dim giving initial vertex location.
//'               Defaults to random locations in a -1,1 square.
//' @param center Numeric vector of length equal to dim specifying the location of the plot center.
//'               Defaults to origin (\code{rep( 0, length( dim ) )}).
//'
//' @return A numeric matrix with as many rows as vertices and as many columns as dim giving the
//'         locations of the vertices in the requested R^dim space.
//'
//' @author José Tomás Atria \email{jtatria@@gmail.com}
//'
//' @family graph layout
//' @export
// [[Rcpp::export]]
RMatD forceatlas(
    RMatD m, int iter=100, bool linlog=false, bool nohubs=false,
    scalar k=1.0, scalar G=1,
    bool strong=false, scalar ks=0.1, scalar ksmax=10.0, scalar delta=0.0, scalar tol=0.1,
    ind dim=2, Nullable<RMatD> init=R_NilValue, Nullable<RVecD> center=R_NilValue
) {
    Mat W     = as_rowmat( m );
    ind n     = W.rows();
    Vec orig  = ( center.isNull() ) ? Vec::Zero( dim ) : as<Vec>( center.get() );
    Mat pos   = ( init.isNull() ) ? Mat::Random( W.rows(), dim ) * 1000: as_rowmat( init.get() );
    Vec deg   = ( W.array() != 0 ).select( Mat::Ones( n, n ), Mat::Zero( n, n ) ).rowwise().sum();
    Mat force = Mat::Zero( n, dim );
    Mat disp  = Mat::Zero( n, dim );
    for( int e = 0; e < iter; e++ ) {
        Mat last = force.eval();
        force = force_func( pos, deg, W, k, delta, G, linlog, nohubs, strong, orig );
        Vec speed = speed_func( force, last, deg, tol, ks, ksmax );
        if( speed.isZero() ) break;
        disp_func( force, speed, disp );
        pos += disp;
    }
    return wrap_rowmat( pos );
}
