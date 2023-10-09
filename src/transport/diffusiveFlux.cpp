#include "diffusiveFlux.h"
#include "lagrange_tmp.h"

/**!
 * Compute values of D(u) along the normal line.
 */
inline std::vector<double> evaluateDiffusionFcns(const vector<double>& uR){

    vector<double> work(samplePoints.size(), 0);

    for (int it = 0; it<work.size(); it++){
        work.at(it) = diffFunc(uR.at(it));   
    }

    return work;
}

inline vector<vertex> assignPts(const vertex& unitNormal,
                                const vertex& mapped,
                                const double& dx){
    int numPts = points.size();

}

double getDiffusiveFluxInterior(const MLWENOUse& mlu,
                                const MeshInfo& mi,
                                const vertexSet& edge,
                                const vertex& unitNormal,
                                const double& len,
                                const indice& globalCellL,
                                const indice& globalCellR,
                                const int& locationL,
                                const int& locationR,
                                const valarray<double>& gwe,
                                const valarray<double>& gpe,
                                const double& scale){

    double work = 0.0;

    int degree = gwe.size() + 1; 
    const int numPts = std::ceil((degree+1)/2.0) * 2;

    // Get diameter
    const double hL = mi.cellArea.at(FlatInidc(mi,globalCellL));
    const double hR = mi.cellArea.at(FlatIndic(mi,globalCellR));
    
    const double h = scale * ((hL < hR) ? hL : hR);

    const double dx = h /(double)(numPts - 1);

    LagrangeBasisDeriv lagDer(numPts - 1);

    const int halfPts = numPts/2;

    // Diffusion function are evaluated at the points in the order of 
    // from outside to inside.
    vector<double> diffVals(numPts, 0);

    for (int g=0; g<gwe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        for (int i=0; i<numPts; i++){
            work -= lagDer.middle(numPts-1, i) / dx * diffVals.at(i) * gwe.at(g) * len/2.0;
        }
    }

    return work;
}
