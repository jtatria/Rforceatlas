#include "work.hpp"
#include "funcs.hpp"

Fa2Worker::Fa2Worker( Mat& pos, const Vec& orig, const GraphData& gd, const Fa2Params& params )
    : pos( pos ), orig( orig ), gd( gd ), params( params )
{
    last    = Mat::Zero( pos.rows(), pos.cols() );
    speed   = 1.0;
    speed_k = 1.0;
}

void Fa2Worker::global_f( Mat& force ) {
    for( ind i = 0; i < pos.rows(); i++ ) {
        Vec vi = pos.row( i );
        // repulsion
        for( ind j = 0; j < i; j++ ) {
            Vec vj = pos.row( j );
            scalar d = dist_func( vi, vj );
            bool overlap = gd.overlap( d, i, j );
            scalar r = repl_func( params.k, d, gd.mass( i ), gd.mass( j ), overlap );
            force.row( i ) += ( vi - vj ) * r;
            force.row( j ) += ( vj - vi ) * r;
        }
        // gravity
        scalar d = dist_func( orig, vi );
        scalar g = grav_func( d, gd.mass( i ), params.G, params.strong );
        force.row( i ) += ( orig - vi ) * g;
    }
}

void Fa2Worker::network_f( Mat& force ) {
    SpMat W = gd.weights();
    for( ind e = 0; e < W.outerSize(); e++ ) {
        for( SpInIt it( W, e ); it; ++it ) {
            ind i = it.row();
            ind j = it.col();
            if( i >= j ) continue;
            Vec vi = pos.row( i );
            Vec vj = pos.row( j );
            scalar d = dist_func( vi, vj );
            scalar w = weight_func( it.value(), params.delta );
            bool overlap = gd.overlap( d, i, j );
            scalar a = attr_func( d, gd.mass( i ), w, overlap, params.linlog, params.nohubs );
            force.row( i ) += ( vi - vj ) * ( a );
            force.row( j ) += ( vj - vi ) * ( a );
        }
    }
}

void Fa2Worker::adjust_speed( const Vec& swing, const Vec& tract ) {
    scalar swing_g = swing.cwiseProduct( ( gd.masses() ).matrix() ).sum();
    scalar tract_g = tract.cwiseProduct( ( gd.masses() ).matrix() ).sum();

    scalar est_jit = 0.05 * std::sqrt<ind>( swing.size() );
    scalar min_jit = std::sqrt<scalar>( est_jit ).real();
    scalar max_jit = 10.0;
    scalar jit = params.tol * std::max( min_jit, std::min( max_jit, est_jit * tract_g / std::pow<scalar>( swing.size(), 2 ) ) );

    scalar min_ks = 0.05;

    if( ( swing_g / tract_g ) > 2 ) {
        if( speed_k > min_ks ) {
            speed_k *= 0.5;
        }
        jit = std::max( jit, params.tol );
    }

    scalar tgt_speed = jit * speed_k * tract_g / swing_g;

    if( swing_g > jit * tract_g ) {
        if( speed_k > min_ks ) speed_k *= 0.7;
    } else if( speed < 1000 ) {
        speed_k *= 1.3;
    }

    speed += std::min( tgt_speed - speed, 0.5 * speed );
}

void Fa2Worker::apply_forces( const Mat& forces, const Vec& swing ) {
    for( ind i = 0; i < pos.rows(); i++ ) {
        scalar factor = speed / ( 1.0 + std::sqrt<scalar>( speed * swing[i] * ( gd.mass( i ) ) ).real() );
        pos.row( i ) += forces.row( i ) * factor;
    }
}

void Fa2Worker::fa2_epoch() {
    Mat force = Mat::Zero( pos.rows(), pos.cols() );
    global_f( force );
    network_f( force );
    Vec swing = swing_vec( force, last );
    Vec tract = tract_vec( force, last );
    adjust_speed( swing, tract );
    apply_forces( force, swing );
    last = force;
}
