#include "diffusiveflux.h"

inline int getsamples(const indice& gcellin,
                      const indice& gcellout,
                      multilevel& ml,
                      mluse& use,
                      double ** lu,
                      vector<double>& samples,
                      const double& dx,
                      const vertex& unitNormal,
                      const vertex& mapped,
                      const Tensor<weights>& allwgts,
                      const std::string& loc){

    // We mark the direction of where unirnormal vector at the edge points to
    // as the out direction.
    // We mark the opposite direction as the in direction

    int halfnumPts = samples.size()/2;

    for (int i=0; i<halfnumPts; i++){

        vertex point0 = mapped - (halfnumPts - 0.5)*dx*unitNormal;
        vertex point1 = mapped + 0.5*dx*unitNormal;

        // Reconstruct values at the given sample points
        double u0 = use.eval(point0, ml, loc, 
                             allwgts({gcellin[0], gcellin[1]}), gcellin, lu);

        double u1 = use.eval(point1, ml, loc, 
                             allwgts({gcellout[0], gcellout[1]}), gcellout, lu);


        samples.at(i)            = diffunc(u0);
        samples.at(i+halfnumPts) = diffunc(u1);

    }

    return 1;
}

double edgefluxintegral(const MeshInfo& mi,
                        const indice& gcellin,
                        const indice& gcellout,
                        const vertexSet& edge,
                        const Tensor<weights>& allwgts,
                        multilevel& ml,
                        mluse& use,
                        double ** lu,
                        const std::string& loc){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double work = 0.0;
    int degree = gwe.size() + 1;

    // number of sample points are designed to be an even number
    const int numPts = std::ceil((degree+1)/2.0) * 2;

    // Compute geometry constant hin and hout
    const double hIn  = mi.cellArea.at(FlatIndic(mi,gcellin));
    const double hOut = mi.cellArea.at(FlatIndic(mi,gcellout));
    
    const double h = 2.0 * ((hIn < hOut) ? hIn : hOut);

    // Compute sample interval
    const double dx = h /(double)(numPts - 1);

    // Compute lagrange interpolation
    LagrangeBasisDeriv lagDer(numPts - 1);

    const int halfPts = numPts/2;

    // Evaluate function at sampling points
    vector<double> samples(numPts, 0);

    // Compute edge related properties
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge, len);

    // Numerical integral on the given edge
    for (int g=0; g<gwe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        getsamples(gcellin, gcellout, ml, use, lu, samples, dx, unitNormal, mapped, allwgts, loc);

        for (int i=0; i<numPts; i++){
            work -= lagDer.middle(numPts-1, i) / dx * gwe[g] * len/2.0 * samples.at(i);
        }        
   
    }
    return work;
}
