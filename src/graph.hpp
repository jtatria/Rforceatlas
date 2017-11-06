#ifndef _GRAPH_HPP
#define _GRAPH_HPP

#include "eigen_types.hpp"

class GraphData {
    const SpMat& W;
    const Vec& size_v;
    Vec mass_v;

public:
    GraphData( const SpMat& W, const Vec& sizes );

    const SpMat& weights() const;

    const Vec& masses() const;

    const Vec& sizes() const;

    const scalar weight( ind i, ind j ) const;

    const scalar mass( ind i ) const;

    const scalar size( ind i ) const;

    const bool overlap( scalar d, ind i, ind j ) const;
};

#endif
