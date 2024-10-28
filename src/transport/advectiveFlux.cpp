#include "advectiveFlux.h"

/**!
 * A Lax-Friedrich style numerical flux scheme
 */
inline double LFFlux(const double& uIn, const double& uOut,
                     const double& fIn, const double& fOut,
                     const vector<double>& coeff){

    return 0.5*(fOut + fIn - coeff[0]*(uOut - uIn));
}

/**!
 * Compute the advective flux at a given point.
 * return a one-sided flux. 
 */
inline double getAdvFluxPoint(const MLWENO::MLWENOUse& mlu, const MeshInfo& mi,
                              const double& u, const vertex& unitNormal){

   std::array<double,2> work = advFunc(u);

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
double getAdvFlux(const MLWENO::MLWENOUse& mlu,
                  const MeshInfo& mi,
                  const std::array<vertex,2>& edge,
                  const vertex& unitNormal,
                  const double& len,
                  const indice& globalCellIn,
                  const indice& globalCellOut,
                  const int& locationIn,
                  const int& locationOut,
                  const valarray<double>& gwe,
                  const valarray<double>& gpe,
                  const vector<double>& coeff){

    std::array<double,2> InFlux;
    std::array<double,2> OutFlux;

    InFlux  = getAdvFluxEdge(mlu, mi, edge, unitNormal, len, globalCellIn , 
                             locationIn , gwe, gpe);
    OutFlux = getAdvFluxEdge(mlu, mi, edge, unitNormal, len, globalCellOut, 
                             locationOut, gwe, gpe);

    return LFFlux(InFlux[0], OutFlux[0], InFlux[1], OutFlux[1], coeff);
}

/**!
 * Compute advective flux on the boundary.
 * Different boundary types require different way of implementation.
 */

double getAdvFlux(const MLWENO::MLWENOUse& mlu,
                  const MeshInfo& mi,
                  const std::array<vertex,2>& edge,
                  const vertex& unitNormal,
                  const double& len,
                  const indice& globalCellIn,
                  const int& locationIn,
                  const valarray<double>& gwe,
                  const valarray<double>& gpe,
                  const vector<double>& coeff,
                  bndryTypeTrans bt){

    // Influx is calculated 
    // Only current targeted element is needed.
    std::array<double,2> InFlux =  getAdvFluxEdge(mlu, mi, edge, unitNormal, len, globalCellIn , 
                                   locationIn , gwe, gpe);

    std::array<double,2> OutFlux;

    double edgeFlux = 0.0;

    switch(bt){
        case wall: 
            OutFlux[0] = -1*InFlux[0];
            OutFlux[1] = -1*InFlux[1];

            edgeFlux = LFFlux(InFlux[0], OutFlux[0], InFlux[1], OutFlux[1], coeff);
        break;

        case freeFlow:

            OutFlux[0] = InFlux[0];
            OutFlux[1] = InFlux[1];

            edgeFlux = LFFlux(InFlux[0], OutFlux[0], InFlux[1], OutFlux[1], coeff);
        break;

        case dirichletAdv:

            OutFlux = bndryValAdv();

            edgeFlux = LFFlux(InFlux[0], OutFlux[0], InFlux[1], OutFlux[1], coeff);
        break;

        case flux:

            edgeFlux = bndryFluxAdv();

        break;

        default : 
            PetscPrintf(PETSC_COMM_WORLD,"Boundary type not prescribed at cell (%d,%d). \n",
                        globalCellIn[0], globalCellIn[1]);
        break;
    }
        return edgeFlux;
}

// =========== Implicit =================================
