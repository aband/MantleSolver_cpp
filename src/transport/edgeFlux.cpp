#include "edgeFlux.h"

double edgeFlux(const valarray<>
){


}

double * edgeFluxAll(const MeshInfo* mi,

                     fluxFunc      fluxfunc,
                     fluxFuncBndry fluxfuncbndry){


    double * edgeFlux = (double *)malloc(sizeof(double) * 
                        mi.MPIlocalHoriEdgeSize * mi.MPIlocalVertEdgeSize);

    // Loop through the entire local mesh chunk
    // Local horizontal edges are looped first
    for (int j=0; j<mi.MPIlocalVertexSize[1]; j++){
        for (int i=0; i<mi.MPIlocalCellSize[0]; i++){

            // Check if it is on the boundary


            // Flaten to integer
            indice local {i,j};
            int flatlocal = FlatIndic(mi.MPIlocalCellSize[0], local);

            // Identify location with global indices 
            indice global = local + mi.MPIlocalCellStart;

            

        }
    }



    return flux;
}
