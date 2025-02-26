#include "diffusiveflux.h"

double getfluxintegral(const MeshInfo& mi,
                       const indice& gcellleft,
                       const indice& gcellright){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double work = 0.0;
    int degree = gwe.size() + 1;
    const int numPts = std::ceil((degree+1)/2.0) * 2;



    return work;
}
