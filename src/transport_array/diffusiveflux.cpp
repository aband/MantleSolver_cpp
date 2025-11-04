#include "diffusiveflux.h"

double edgefluxintegral(const reconstruction& recon_neg,
                        const reconstruction& recon_pos,
								const indice& gcellin,
								const indice& gcellout,
                        const vector<tensorstencilpoly>& sten_lg,
                        const vector<tensorstencilpoly>& sten_sm,
                        const int sampleSize, 
                        double ** localvals,
                        const vertexSet& edge){

    double work = 0.0;

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);
    vertex unitNormal = UnitNormal(edge, len);

    int degree = gwe.size() + 1;

    const int halfPts = std::ceil((degree+1)/2.0);
    const int numPts =  halfPts * 2;

    const double hIn  = mi.cellArea.at(FlatIndic(mi,gcellin));
    const double hOut = mi.cellArea.at(FlatIndic(mi,gcellout));

    const double h = 2.0 * ((hIn < hOut) ? hIn : hOut);

    // Compute sample interval
    const double dx = h /(double)(numPts - 1);

    // Compute lagrange interpolation
    LagrangeBasisDeriv lagDer(numPts - 1);

    for (int g=0; g<gwe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g], edge});



		  for (int i=0; i<numPts; i++){
            work -= lagDer.middle(numPts-1, i)/dx * gwe[g] * len/2.0 * samples.at(i);
		  } 

    }

    return work;
}
