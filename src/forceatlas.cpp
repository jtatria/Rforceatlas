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

//' @backref src/forceatlas.cpp
#include "Rforceatlas_types.hpp"

using namespace Rcpp;

Mat dist_func( Mat pos, scalar min=0.01 ) {
    Mat dist( pos.rows(), pos.rows() );
    for( ind i = 0; i < pos.rows(); i++ ) {
        for( ind j = 0; j < pos.rows(); j++ ) {
            Vec vi = pos.row( i );
            Vec vj = pos.row( j );
            dist.coeffRef( i, j ) = ( vi - vj ).norm();
        }
    }
    dist = ( min > 0 ) ? ( dist.array() < min ).select( min, dist ) : dist;
    return dist;
}

Mat repl_func( Vec deg, Mat dist, scalar k=400 ) {
    Mat out = ( ( deg.array() + 1 ).matrix() * ( deg.array() + 1 ).matrix().transpose() )
        .cwiseQuotient( dist ) * k;
    return out.eval();
}

Mat attr_func(
    Mat dist, Mat wgts, Vec deg, scalar delta=1.0, bool linlog=false, bool nohubs=false
) {
    Mat out = linlog ? dist.array().log1p() : dist;
    out = delta != 0 ? out.cwiseProduct( wgts.pow( delta ) ) : out;
    out = nohubs ? out.cwiseProduct( deg * Vec::Ones( deg.size() ).transpose() ) : out;
    return out;
}

Mat grav_func( Mat pos, Vec deg, Vec orig=Vec::Zero( 2 ), scalar G=1.0 ) {
    Mat out( pos.rows(), pos.cols() );
    for( ind i = 0; i < pos.rows(); i++ ) {
        Vec uv = ( orig.transpose() - pos.row( i ) );
        uv = ( uv.array() / uv.norm() * ( G * ( deg( i ) + 1.0 ) ) );
        out.row( i ) = uv;
    }
    return out;
}

Mat force_func( Mat attr, Mat repl, Mat grav, Mat pos, Mat dist ) {
    Mat out( pos.rows(), pos.cols() );
    for( ind i = 0; i < pos.cols(); i++ ) {
        Mat m1 = pos.col( i ) * Vec::Ones( pos.rows() ).transpose();
        Mat m2 = Vec::Ones( pos.rows() ) * pos.col( i ).transpose();
        out.col( i ) = ( ( m1 - m2 ).cwiseQuotient( dist ).cwiseProduct( repl - attr ) )
            .rowwise().sum();
    }
    out = out + grav;
    return out;
}

Vec speed_func( Mat f_t0, Mat f_t1, Vec deg, scalar tol=.1, scalar s=.1, scalar smax=10.0 ) {
    Vec swing = ( f_t0 - f_t1 ).rowwise().sum().cwiseAbs();
    if( ( swing.array() == 0.0 ).all() ) return swing;
    Vec tract = ( f_t0 + f_t1 ).rowwise().sum().cwiseAbs();
    scalar swing_g = swing.cwiseProduct( ( deg.array() + 1.0 ).matrix() ).sum();
    scalar tract_g = tract.cwiseProduct( ( deg.array() + 1.0 ).matrix() ).sum();
    scalar speed_g = tol * tract_g / swing_g;
    Vec speed_max = ( f_t0.rowwise().norm() ).array().pow( -1.0 ) * smax;
    Vec speed = ( swing.cwiseSqrt() * ( speed_g + 1.0 ) ).array().pow( -1 ) * ( s * speed_g );
    speed = ( speed.array() > speed_max.array() ).select( speed_max, speed );
    return speed;
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
//'               Defaults to random locations in a -1000,1000 square.
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
    RMatD m, int iter=100, bool linlog=false, bool nohubs=false, scalar k=400, scalar G=1,
    scalar ks=0.1, scalar ksmax=10.0, scalar delta=1.0, scalar tol=0.1, ind dim=2,
    Nullable<RMatD> init=R_NilValue, Nullable<RVecD> center=R_NilValue
) {
    Mat A      = as<Mat>( m );
    ind n      = A.rows();
    Vec orig   = ( center.isNull() ) ? Vec::Zero( dim ) : as<Vec>( center.get() );
    Mat pos    = ( init.isNull() ) ? Mat::Random( A.rows(), dim ) * 1000 : as<Mat>( init.get() );
    Vec deg    = ( A.array() != 0 ).select( Mat::Ones( n, n ), Mat::Zero( n, n ) ).rowwise().sum();
    Mat force = Mat::Zero( n, dim );
    for( int e = 0; e < iter; e++ ) {
        Mat last = force.eval();
        Mat dist = dist_func( pos );
        Mat repl = repl_func( deg, dist, k );
        Mat attr = attr_func( dist, A, deg, delta, linlog, nohubs );
        Mat grav = grav_func( pos, deg, orig, G );
        force = force_func( attr, repl, grav, pos, dist );
        Vec speed = speed_func( force, last, deg, tol, ks, ksmax );
        if( ( speed.array() == 0.0 ).all() ) break;
        for( ind i = 0; i < force.rows(); i++ ) {
            pos.row( i ) += force.row( i ) * speed( i );
        }
    }
    return wrap( pos );
}
