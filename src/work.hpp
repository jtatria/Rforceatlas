#ifndef _WORK_HPP
#define _WORK_HPP 1

#include "eigen_types.hpp"
#include "graph.hpp"
#include "params.hpp"

class Fa2Worker {
    Mat& pos;

    const Vec& orig;
    const GraphData& gd;
    const Fa2Params params;

    scalar speed;
    scalar speed_k;
    Mat last;

    void global_f( Mat& force );
    void network_f( Mat& force );
    void adjust_speed( const Vec& swing, const Vec& tract );
    void apply_forces( const Mat& forces, const Vec& swing );

public:
    Fa2Worker( Mat& pos, const Vec& orig, const GraphData& gd, const Fa2Params& params );
    void fa2_epoch();
};

#endif
