#include "funcs.hpp"

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
