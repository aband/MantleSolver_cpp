#include "edgeFlux.h"

/**!
 * Update edge flux on a given edge.
 * cell "globalCellIn" and cell "globalCellOut" sharing this edge.
 */
inline double getEdgeFluxBoundary(const MLWENOUse& mluAdv,
                                  const MLWENOUse& mluDif,
                                  const MeshInfo& mi,
                                  const vertexSet& edge,
                                  const vertex& unitNormal,
                                  const double& len,
                                  const indice& globalCellIn,
                                  const indice& globalCellOut,
                                  const int& locationIn,
                                  const int& locationOut,
                                  const valarray<double>& gwe,
                                  const valarray<double>& gpe,
                                  const double& alpha,
                                  const double& scale,
                                  const std::string& boundaryTypeAdv,
                                  const std::string& boundaryTypeDif,){

    double work = 0.0;

    return work;
}

void Transport::edgeFlux(const MeshInfo& mi,
                         const std::array& flowType){

    edgeFlux_.resize(mi.MPIlocalHoriEdgeSize + 
                     mi.MPIlocalVertEdgeSize);

    // Assign flow type to the class edgeflux
    switch(flowType){
        case "advection":
            isAdv = true;
            break;
        case "diffusion":
            isDif = true;
            break;
        case "adv-diff":
            isAdv = true;
            isDif = true;
            break;
        default:
            cout << "Flow type not declared !" << endl;
            break;
    }
}

/**!
 * Update flux on all the edges.
 */
void Transport::getEdgeFlux(){
    // Initialize map holding all edge flux.
}
