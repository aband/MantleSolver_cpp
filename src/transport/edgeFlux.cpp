#include "edgeFlux.h"

inline bool outsideBoundary(const MeshInfo& mi,
                            const indice& cell){
    if (cell[0]<0 || cell[0]>mi.MPIglobalCellSize[0] || 
        cell[1]<0 || cell[1]>mi.MPIglobalCellSize[1]){
        return 0;
    } else {
        return true;
    }

}

// Check if the edge is on boundary or not.
inline bool onBoundary(const MeshInfo& mi,
                       const std::array<indice, 2>& nbr){

    // Check if any cell is outside of the boundary
    if (outsideBoundary(mi,nbr.at(0)) || outsideBoundary(mi,nbr.at(1))){
        return true;
    } else {
        return false;
    }
}

EdgeFlux::EdgeFlux(const MeshInfo& mi,
                   const flowType& fT){

    edgeFlux_.resize(mi.MPIlocalHoriEdgeSize + 
                     mi.MPIlocalVertEdgeSize);

    // Assign flow type to the class edgeflux
    switch(fT){
        case advection:
            isAdv = true;
            break;
        case diffusion:
            isDif = true;
            break;
        case adv_diff:
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
 * Distinguish between different boundary condition here.
 */
void EdgeFlux::getEdgeFlux(const MeshInfo& mi,
                           const MLWENO::MLWENOUse& mluAdv,
                           const MLWENO::MLWENOUse& mluDif){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Simplified example, where only interior location is identified
    int locationL = 0;
    int locationR = 0;

    // Compute edge flux
    for (int e=0 ;e<edgeFlux_.size(); e++){
        // Switch local edge index to global edge index.
        const int globalEdge = edgeIndexLocalToGlobal(mi,e); 

        // Extract two cell index sharing the given edge.
        // Cells are given in global cell indice.
		  const std::array<indice,2> nBrs = extractEdgeNbr(mi, globalEdge);

        // Assign flux values to edges
        edgeFlux_.at(e) = 0;

        const std::array<vertex, 2> edge = extractEdge(mi, globalEdge);

        double len = getEdgeLength(edge); 

        vertex unitNormal = getUnitNormal(edge,len);

        // Judging whether the edge is on the boundary or not.
        if (onBoundary(mi,nBrs)){
            // Check if inflow or outflow
				// No flow condition
            edgeFlux_.at(e) = 0; 

        } else { // Interior

            if (isAdv == 1){
                edgeFlux_.at(e) += getAdvFluxInterior(mluAdv, mi, 
                                     edge, unitNormal, len, nBrs[1], nBrs[0], 
                                     locationL, locationR, gwe, gpe, alpha_);
            }

            if (isDif == 1){
                edgeFlux_.at(e) += getDifFluxInterior(mluDif, mi, 
                                     edge, unitNormal, len, nBrs[1], nBrs[0], 
                                     locationL, locationR, gwe, gpe, scale_);
            }

        }
    }
}
