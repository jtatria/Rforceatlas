#ifndef _FA2_HPP
#define _FA2_HPP 1

#include "eigen_types.hpp"

struct Fa2Params {
public:
    const scalar delta;
    const scalar tol;
    const scalar k;
    const scalar G;
    const bool   linlog;
    const bool   strong;
    const bool   nohubs;
    const bool   overlap;

    Fa2Conf(
        const scalar delta=1.0,
        const scalar tol=1.0,
        const scalar k=10.0,
        const scalar G=1.0,
        const bool linlog=false,
        const bool strong=false,
        const bool nohubs=false,
        const bool overlap=false
    ) :
        delta( delta ),
        tol( tol ),
        k( k ),
        G( G ),
        linlog( linlog ),
        strong( strong ),
        nohubs( nohubs ),
        overlap( overlap )
    {}
};

#endif
