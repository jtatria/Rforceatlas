#include <Rcpp.h>
#include "Rforceatlas_types.h"
#include "params.hpp"
#include "graph.hpp"
#include "work.hpp"

using namespace Rcpp;

inline Mat as_rowmat( RMatD input ) {
    Eigen::MatrixXd tmp = Rcpp::as<Eigen::MatrixXd>( input );
    Mat out( tmp.rows(), tmp.cols() );
    out = tmp;
    return out;
}

inline RMatD wrap_rowmat( Mat output ) {
    Eigen::MatrixXd tmp( output.rows(), output.cols() );
    tmp = output;
    return Rcpp::wrap( tmp );
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
    S4 m, Nullable<RMatD> init=R_NilValue, Nullable<RVecD> center=R_NilValue, Nullable<RVecD> vsizes=R_NilValue,
    int dim=2, int iter=100, scalar delta=1.0, scalar tol=1.0, scalar k=10.0, scalar G=1.0,
    bool linlog=false, bool strong=false, bool nohubs=false, bool overlap=false
) {
    Fa2Params params( delta, tol, k, G, linlog, strong, nohubs, overlap );
    SpMat W   = as<SpMat>( m );
    Vec orig  = ( center.isNull() ) ? Vec::Zero( dim ) : as<Vec>( center.get() );
    Mat pos   = ( init.isNull() ) ? Mat::Random( W.rows(), 2 ) * 1000: as_rowmat( init.get() );
    Vec sizes  = ( vsizes.isNull() ) ? Vec::Ones( W.rows() ) : as<Vec>( vsizes );
    GraphData gd( W, sizes );
    Fa2Worker wrkr( pos, orig, gd, params );
    for( int i = 0; i < iter; i++ ) {
        wrkr.fa2_epoch();
    }
    return wrap_rowmat( pos );
}
