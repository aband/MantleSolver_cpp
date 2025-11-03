#include "diffusiveflux.h"

double edgefluxintegral(const reconstruction& recon_neg,
                        const reconstruction& recon_pos,
                        const vector<tensorstencilpoly>& sten_lg,
                        const vector<tensorstencilpoly>& sten_sm,
                        double ** localvals,
                        const vertexSet& edge){

    double work = 0.0;

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);
    vertex unitNormal = UnitNormal(edge, len);



    return work;
}
