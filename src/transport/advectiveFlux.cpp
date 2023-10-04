#include "advectionFlux.h"

/**!
 * A Lax-Friedrich style numerical flux scheme
 */
inline double numericalFlux(const double& uIn, const double& uOut,
                            const double& fIn, const double& fOut,
                            const double& alpha){

    return 0.5*(fIn + fOut - alpha*(uOut - uIn));
}

/**!
 * Compute the advective flux at a given point.
 * return a one-sided flux. 
 */
inline double getAdvFluxPoint(const MLWENOUse& mlu, const MeshInfo& mi,
                              const double& uR, const vertex& unitNormal){

   vertex work {funcX(uR), 
                funcY(uR)}; 

   return std::inner_product(work.begin(), work.end(), unitNormal.begin(), 0);
}

/**!
 * Integrate one sided flux along the edge.
 */
inline std::array<double,2> getAdvFluxEdge(const MLWENOUse& mlu, 
                                           const MeshInfo& mi, 
                                           const vertexSet& edge, 
                                           const vertex& unitNormal,
                                           const double& len,
                                           const indice& globalCell,
                                           const int& location,
                                           const valarray<double>& gwe,
                                           const valarray<double>& gpe){
    double work1 = 0.0;
    double work2 = 0.0;

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]},edge);
        // Get reconstructed value at the given point
        uR = mlu.Evaluate(point,globalCell,mi,location);
        work1 += gwe[g]*getAdvFluxPoint(mlu,mi,uR,unitNormal) * len/2.0;
        work2 += gwe[g]*uR * len/2.0;
    }

    return {work1, work2};
}

/**!
 * Compute advective interior of the domain.
 */

inline double getAdvFluxInterior(const MLWENOUse& mlu,
                                 const MeshInfo& mi,
                                 const vertexSet& edge,
                                 const vertex& unitNormal,
                                 const double& len,
                                 const indice& globalCellIn,
                                 const indice& globalCellOut,
                                 const int& location,
                                 const valarray<double>& gwe,
                                 const valarray<double>& gpe,
                                 const double& alpha){

    std::array<double,2> In;
    std::array<double,2> Out;

    In  = getAdvFluxEdge(mlu, mi, edge, unitNormal, len, globalCellIn , location, gwe, gpe);
    Out = getAdvFluxEdge(mlu, mi, edge, unitNormal, len, globalCellOut, location, gwe, gpe);

    return numericalFlux(In[0], Out[0], In[1], Out[1], alpha);
}

/**!
 * Compute advective flux on the boundary.
 * Inflow and outflow boundary are discussed separately.
 */
inline double getAdvFluxBoundary(const std::string& inflowType,
                                 const std::string& outflowType,
                                 const vertex& unitNormal){

    switch
}

double getAdvFlux(const MLWENOUse& mlu,
                  const MeshInfo& mi,
                  const std::string& ){

}
