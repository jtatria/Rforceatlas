#include "eigen_types.hpp"

bool get_overlap( const scalar d, const ind i, const ind j ) {
    return false;
}

scalar edge_weight( scalar delta, scalar w ) {
    if( delta <= 0.0 ) return 1.0;
    if( delta <= 1.0 ) return w;
    return std::pow<scalar>( w, delta );
}

Vec deg_vec( const SpMat& W ) {
    Vec deg( W.outerSize() );
    for( ind i = 0; i < W.outerSize(); i++ ) {
        deg[i] = W.row( i ).sum();
    }
    return deg;
}

scalar attr_func(
    const scalar d, const scalar deg_i, const scalar wgt,
    const bool linlog, const bool nohubs, const bool overlap
) {
    if( overlap ) return 0;
    scalar mass = deg_i + 1.0;
    scalar a = wgt;
    a *= linlog ? std::log<scalar>( 1.0 + d ).real() / d : 1.0;
    a /= nohubs ? mass : 1.0;
    return a;
}

scalar repl_func(
    const scalar k, const scalar d, const scalar deg_i, const scalar deg_j,
    const bool overlap
) {
    scalar mass = ( deg_i + 1.0 ) * ( deg_j + 1.0 );
    scalar r = k * mass * ( overlap ? 100 : 1 / d / d );
    return r;
}

scalar grav_func(
    const scalar d, const scalar deg_i, const scalar G, const bool strong
) {
    scalar mass = ( deg_i + 1.0 );
    scalar g = ( G * mass ) / d;
    return g;
}

scalar dist_func( const Vec& vi, const Vec& vj ) {
    return ( vi - vj ).norm();
}

void global_forces( Mat& force, const Mat& pos, const Vec& deg, const scalar k, const scalar G, const Vec& orig, const bool strong ) {
    for( ind i = 0; i < pos.rows(); i++ ) {
        // gravity
        Vec vi = pos.row( i );
        scalar d = dist_func( orig, vi );
        scalar g = grav_func( d, deg[i], G, strong );
        force.row( i ) = ( orig - vi ) *g;
        // repulsion
        for( ind j = 0; j < pos.rows(); j++ ) {
            if( i == j ) continue;
            Vec vj = pos.row( j );
            scalar d = dist_func( vi, vj );
            bool overlap = get_overlap( d, i, j );
            scalar r = repl_func( k, d, deg[i], deg[j], overlap );
            force.row( i ) += ( vi - vj ) * r;
        }
    }
}

void network_forces(
        Mat& force, const Mat& pos, const SpMat& W, const Vec& deg,
        const scalar delta, const bool linlog, const bool nohubs
) {
    for( ind e = 0; e < W.outerSize(); e++ ) {
        for( SpInIt it( W, e ); it; ++it ) {
            ind i = it.row();
            ind j = it.col();
            if( i == j ) continue;
            Vec vi = pos.row( i );
            Vec vj = pos.row( j );
            scalar d = dist_func( vi, vj );
            scalar w = edge_weight( delta, it.value() );
            bool overlap = get_overlap( d, i, j );
            scalar a = attr_func( d, deg[i], w, linlog, nohubs, overlap );
            force.row( i ) += ( vi - vj ) * ( a * -1.0 );
        }
    }
}

void displacement(
        Mat& disp, const Mat& f_t0, const Mat& f_t1, const Vec& deg,
        const scalar tol, scalar& speed, scalar& speed_k
) {
    ind n = f_t0.rows();

    Vec swing = ( f_t1 - f_t0 ).rowwise().norm();
    Vec tract = ( f_t1 + f_t0 ).rowwise().norm() / 2;
    scalar swing_g = swing.cwiseProduct( ( deg.array() + 1.0 ).matrix() ).sum();
    scalar tract_g = tract.cwiseProduct( ( deg.array() + 1.0 ).matrix() ).sum();

    scalar est_jit = 0.05 * std::sqrt<ind>( n );
    scalar min_jit = std::sqrt<scalar>( est_jit ).real();
    scalar max_jit = 10.0;
    scalar jit = tol * std::max( min_jit, std::min( max_jit, est_jit * tract_g / std::pow<ind>( n, 2 ) ) );

    scalar min_ks = 0.05;
    if( ( swing_g / tract_g ) > 2 ) {
        if( speed_k > min_ks ) {
            // GLOBAL!
            speed_k *= 0.5;
        }
        jit = std::max( jit, tol );
    }

    scalar tgt_speed = jit * speed_k * tract_g / swing_g;
    if( swing_g > jit * tract_g ) {
        if( speed_k > min_ks ) speed_k *= 0.7;
    } else if( speed < 1000 ) {
        speed_k *= 1.3;
    }

    // GLOBAL!
    speed = speed + std::min( tgt_speed - speed, 0.5 * speed );

    for( ind i = 0; i < n; i++ ) {
        scalar swing = deg[i] * ( f_t1.row( i ) - f_t0.row( i ) ).norm();
        scalar factor = speed / ( 1.0 + std::sqrt<scalar>( speed * swing ).real() );
        disp.row( i ) = f_t0.row( i ) * factor;
    }
}

Mat fa2(
    Mat& pos,
    const SpMat& W, const Vec& orig,
    const int iter=100, const ind dim=2,
    const scalar delta=1.0, const scalar tol=1.0, const scalar k=10.0, const scalar G=1.0,
    const bool linlog=false, const bool strong=false, const bool nohubs=false, const bool overlap=false
) {
    Vec deg = deg_vec( W );
    scalar speed = 1.0;
    scalar speed_k = 1.0;
    Mat force = Mat::Zero( W.rows(), dim );
    Mat disp  = Mat::Zero( W.rows(), dim );
    for( int e = 0; e < iter; e++ ) {
        Mat last = force.eval();
        global_forces( force, pos, deg, k, G, orig, strong );
        network_forces( force, pos, W, deg, delta, linlog, nohubs );
        displacement( disp, force, last, deg, tol, speed, speed_k );
        pos += disp;
    }
    return pos;
}
