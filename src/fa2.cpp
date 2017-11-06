#include "eigen_types.hpp"

bool get_overlap( const scalar d, const ind i, const ind j ) {
    return false;
}

scalar edge_weight( scalar delta, scalar w ) {
    if( delta <= 0.0 ) return 1.0;
    if( delta <= 1.0 ) return w;
    return std::pow<scalar>( w, delta );
}

Vec mass_vec( const SpMat& W ) {
    Vec mass = Vec::Zero( W.outerSize() );
    for( ind i = 0; i < W.outerSize(); i++ ) {
        for( SpInIt it( W, i ); it; ++it ) {
            mass[i] += it.value() > 0 ? 1 : 0;
        }
        mass[i] += 1.0;
    }
    return mass;
}

scalar attr_func(
    const scalar d, const scalar mass_i, const scalar wgt,
    const bool linlog, const bool nohubs, const bool overlap
) {
    if( overlap ) return 0;
    scalar mass = mass_i;
    scalar a = wgt;
    a *= linlog ? std::log<scalar>( 1.0 + d ).real() / d : 1.0;
    a /= nohubs ? mass : 1.0;
    return a * -1.0;
}

scalar repl_func(
    const scalar k, const scalar d, const scalar mass_i, const scalar mass_j,
    const bool overlap
) {
    scalar mass = ( mass_i ) * ( mass_j );
    scalar r = k * mass * ( overlap ? 100 : 1 / d / d );
    return r;
}

scalar grav_func(
    const scalar d, const scalar mass_i, const scalar G, const bool strong
) {
    scalar mass = ( mass_i );
    scalar g = ( G * mass ) / ( strong ? 1.0 : d );
    return g;
}

scalar dist_func( const Vec& vi, const Vec& vj ) {
    return ( vi - vj ).norm();
}

void global_forces( Mat& force, const Mat& pos, const Vec& mass, const scalar k, const scalar G, const Vec& orig, const bool strong ) {
    for( ind i = 0; i < pos.rows(); i++ ) {
        Vec vi = pos.row( i );
        // repulsion
        for( ind j = 0; j < i; j++ ) {
            Vec vj = pos.row( j );
            scalar d = dist_func( vi, vj );
            bool overlap = get_overlap( d, i, j );
            scalar r = repl_func( k, d, mass[i], mass[j], overlap );
            force.row( i ) += ( vi - vj ) * r;
            force.row( j ) += ( vj - vi ) * r;
        }
        // gravity
        scalar d = dist_func( orig, vi );
        scalar g = grav_func( d, mass[i], G, strong );
        force.row( i ) += ( orig - vi ) *g;
    }
}

void network_forces(
        Mat& force, const Mat& pos, const SpMat& W, const Vec& mass,
        const scalar delta, const bool linlog, const bool nohubs
) {
    for( ind e = 0; e < W.outerSize(); e++ ) {
        for( SpInIt it( W, e ); it; ++it ) {
            ind i = it.row();
            ind j = it.col();
            if( i >= j ) continue;
            Vec vi = pos.row( i );
            Vec vj = pos.row( j );
            scalar d = dist_func( vi, vj );
            scalar w = edge_weight( delta, it.value() );
            bool overlap = get_overlap( d, i, j );
            scalar a = attr_func( d, mass[i], w, linlog, nohubs, overlap );
            force.row( i ) += ( vi - vj ) * ( a );
            force.row( j ) += ( vj - vi ) * ( a );
        }
    }
}

void displacement(
        Mat& disp, const Mat& f_t0, const Mat& f_t1, const Vec& mass,
        const scalar tol, scalar& speed, scalar& speed_k
) {
    ind n = f_t0.rows();

    Vec swing = ( f_t1 - f_t0 ).rowwise().norm();
    Vec tract = ( f_t1 + f_t0 ).rowwise().norm() / 2;
    scalar swing_g = swing.cwiseProduct( ( mass ).matrix() ).sum();
    scalar tract_g = tract.cwiseProduct( ( mass ).matrix() ).sum();

    scalar est_jit = 0.05 * std::sqrt<ind>( n );
    scalar min_jit = std::sqrt<scalar>( est_jit ).real();
    scalar max_jit = 10.0;
    scalar jit = tol * std::max( min_jit, std::min( max_jit, est_jit * tract_g / std::pow<scalar>( n, 2 ) ) );

    scalar min_ks = 0.05;

    if( ( swing_g / tract_g ) > 2 ) {
        if( speed_k > min_ks ) {
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

    scalar speed_adj = std::min( tgt_speed - speed, 0.5 * speed );
    speed += speed_adj;

    Vec fvec( n );
    for( ind i = 0; i < n; i++ ) {
        scalar factor = speed / ( 1.0 + std::sqrt<scalar>( speed * swing[i] * ( mass[i] ) ).real() );
        fvec[i] = factor;
        disp.row( i ) = f_t0.row( i ) * factor;
    }
}

void apply( Mat& pos, Mat& disp ) {
    pos += disp;
}

Mat fa2(
    Mat& pos,
    const SpMat& W, const Vec& orig,
    const int iter=100, const ind dim=2,
    const scalar delta=1.0, const scalar tol=1.0, const scalar k=10.0, const scalar G=1.0,
    const bool linlog=false, const bool strong=false, const bool nohubs=false, const bool overlap=false
) {
    Vec deg = mass_vec( W );
    scalar speed = 1.0;
    scalar speed_k = 1.0;
    Mat force = Mat::Zero( W.rows(), dim );
    Mat disp  = Mat::Zero( W.rows(), dim );
    for( int e = 0; e < iter; e++ ) {
        Mat last = force;
        force = Mat::Zero( W.rows(), dim );
        global_forces( force, pos, deg, k, G, orig, strong );
        network_forces( force, pos, W, deg, delta, linlog, nohubs );
        displacement( disp, force, last, deg, tol, speed, speed_k );
        apply( pos, disp );
    }
    return pos;
}
