#include "diffusiveFlux.h"

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

double getDiffusiveFluxInterior(const MLWENOUse& mlu,
                                const MeshInfo& mi,
                                const vertexSet& edge,
                                const vertex& unitNormal,
                                const double& len,
                                const indice& globalCellIn,
                                const indice& globalCellOut,
                                const int& location,
                                const valarray<double>& gwe,
                                const valarray<double>& gpe,
                                const double& scale){

    double work = 0.0;

    int degree = gwe.size() + 1; 
    const int numPts = std::ceil((degree+1)/2.0) * 2;

    // Get diameter
    const double hIn = ...
    const double hOut = ...
    
    const double h = scale * ()

    const double dx = h /(double)(numPts - 1);

    for (int g=0; g<gwe.size(); g++){
        vector<double> = 

    }

    return work;
}
