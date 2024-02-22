#include "advectiveFlux.h"

/**!
 * A Lax-Friedrich style numerical flux scheme
 */
inline double numericalFlux(const double& uL, const double& uR,
                            const double& fL, const double& fR,
                            const double& alpha){

    return 0.5*(fR + fL - alpha*(uR - uL));
}

/**!
 * Compute the advective flux at a given point.
 * return a one-sided flux. 
 */
inline double getAdvFluxPoint(const MLWENO::MLWENOUse& mlu, const MeshInfo& mi,
                              const double& uR, const vertex& unitNormal){

   std::array<double,2> work = advFunc(uR);

   return work[0]*unitNormal[0]+work[1]*unitNormal[1];
}

/**!
 * Integrate one sided flux along the edge.
 */
inline std::array<double,2> getAdvFluxEdge(const MLWENO::MLWENOUse& mlu, 
                                           const MeshInfo& mi, 
                                           const std::array<vertex,2>& edge, 
                                           const vertex& unitNormal,
                                           const double& len,
                                           const indice& globalCell,
                                           const int& location,
                                           const valarray<double>& gwe,
                                           const valarray<double>& gpe){
    double work1 = 0.0;
    double work2 = 0.0;

    vertexSet tmpEdge = {edge[0], edge[1]};

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]},tmpEdge);
        // Get reconstructed value at the given point
        double uR = mlu.Evaluate(mapped,globalCell,mi,location);
        work1 += gwe[g]*getAdvFluxPoint(mlu,mi,uR,unitNormal) * len/2.0;
        work2 += gwe[g]*uR * len/2.0;
    }

    return {work1, work2};
}

/**!
 * Compute advective interior of the domain.
 */
double getAdvFluxInterior(const MLWENO::MLWENOUse& mlu,
                          const MeshInfo& mi,
                          const std::array<vertex,2>& edge,
                          const vertex& unitNormal,
                          const double& len,
                          const indice& globalCellL,
                          const indice& globalCellR,
                          const int& locationL,
                          const int& locationR,
                          const valarray<double>& gwe,
                          const valarray<double>& gpe,
                          const double& alpha){

    std::array<double,2> LFlux;
    std::array<double,2> RFlux;

    LFlux = getAdvFluxEdge(mlu, mi, edge, unitNormal, len, globalCellL, locationL, gwe, gpe);
    RFlux = getAdvFluxEdge(mlu, mi, edge, unitNormal, len, globalCellR, locationR, gwe, gpe);

    return numericalFlux(LFlux[0], RFlux[0], LFlux[1], RFlux[1], alpha);
}

// =========== Implicit =================================
