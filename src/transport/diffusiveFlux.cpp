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

inline double getAdvFluxEdge(){

    double work;

    return work;
}

double getDiffusiveFluxInterior(const MLWENOUse& mlu){

    double work = 0.0;

    return work;
}
