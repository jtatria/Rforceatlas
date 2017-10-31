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
    Mat out = ( ( deg.array() + 1 ).matrix() * ( deg.array() + 1 ).matrix().transpose() ).cwiseQuotient( dist ) * k;
    return out.eval();
}

Mat attr_func( Mat dist, Mat wgts, Vec deg, scalar delta=1.0, bool linlog=false, bool nohubs=false ) {
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
        out.col( i ) = ( ( m1 - m2 ).cwiseQuotient( dist ).cwiseProduct( repl - attr ) ).rowwise().sum();
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

RMatD forceatlas2(
    RMatD m, int iter=100, bool linlog=false, bool nohubs=false, scalar k=400, scalar G=1,
    scalar ks=0.1, scalar ksmax=10.0, scalar delta=1.0, scalar tol=0.1, ind dim=2,
    Nullable<RMatD> init_pos=R_NilValue,
    Nullable<RVecD> orig_pos=R_NilValue, int plotstep=10
) {
    Mat A      = as<Mat>( m );
    ind n      = A.rows();
    Vec orig   = ( orig_pos.isNull() ) ? Vec::Zero( dim ) : as<Vec>( orig_pos.get() );
    Mat pos    = ( init_pos.isNull() ) ? Mat::Random( A.rows(), dim ) * 1000 : as<Mat>( init_pos.get() );
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
