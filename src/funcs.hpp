#ifndef _FUNCS_HPP
#define _FUNCS_HPP

#include "eigen_types.hpp"

scalar attr_func(
    const scalar d, const scalar mass_i, const scalar wgt,
    const bool linlog, const bool nohubs, const bool overlap
);

scalar repl_func(
    const scalar k, const scalar d, const scalar mass_i, const scalar mass_j,
    const bool overlap
);

scalar grav_func(
    const scalar d, const scalar mass_i, const scalar G, const bool strong
);

scalar dist_func( const Vec& vi, const Vec& vj );

#endif
