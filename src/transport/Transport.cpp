#include "edgeFlux.h"

inline bool outsideBoundary(const indice& cell){


}

// Check if the edge is on boundary or not.
inline bool onBoundary(const std::array<indice, 2>& nbr){

    if (nbr.at(0)){

        return true;
    }
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
 * Distinguish between different boundary condition here.
 */
void Transport::getEdgeFlux(const MeshInfo& mi){
    // Compute edge flux
    for (int e=0 ;e<edgeFlux_.size(); i++){
        // Switch local edge index to global edge index.
        int globalEdge = edgeIndexLocalToGlobal(mi,e); 

        // Extract two cell index sharing the given edge.
        // Cells are given in global cell indice.
        vertex nBrs = extractEdgeNbr(mi, globalEdge);

        // Assign flux values to edges
        edgeFlux_.at(e) = 0;

        std::array<vertex, 2> edge = extractEdge(mi, globalEdge);

        double len = getEdgeLength(edge); 

        vertex unitNormal = getUnitNormal(edge,len);

        // Judging whether the edge is on the boundary or not.
        if (onBoundary(nBrs)){


        } else {

            if (isAdv == 1){
                edgeFlux_.at(e) += getAdvFluxInterior(mluAdv, mi, globalEdge);
            }

            if (isDif == 1){
                edgeFlux_.at(e) += getDifFluxInterior();
            }

        }
    }
}
