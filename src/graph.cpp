#include "graph.hpp"

GraphData::GraphData( const SpMat& W, const Vec& sizes ) : W( W ), size_v( sizes ) {
    mass_v = Vec::Zero( W.outerSize() );
    for( ind i = 0; i < W.outerSize(); i++ ) {
        for( SpInIt it( W, i ); it; ++it ) {
            mass_v[i] += it.value() > 0 ? 1 : 0;
        }
        mass_v[i] += 1.0;
    }
}

const SpMat& GraphData::weights() const {
    return W;
}

const Vec& GraphData::masses() const {
    return mass_v;
}

const Vec& GraphData::sizes() const {
    return size_v;
}

const scalar GraphData::weight( ind i, ind j ) const {
    return W.coeff( i, j );
}

const scalar GraphData::mass( ind i ) const {
    return mass_v[i];
}

const scalar GraphData::size( ind i ) const {
    return size_v[i];
}

const bool GraphData::overlap( scalar d, ind i, ind j ) const {
    return false;
}
