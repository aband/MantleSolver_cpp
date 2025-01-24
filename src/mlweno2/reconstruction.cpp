#include "reconstruction.h"

int reconstruction::prepare(const vector<int>& insize, const MeshInfo& mi){

    size = insize;

    // Counting stencils by its left bottom cell 
    for (int j=mi.MPIlocalCellStart[1]-mi.cellGhostLayerSize; 
         j<mi.MPIlocalCellStart[1]+mi.MPIlocalCellSize[1]+mi.cellGhostLayerSize; j++){
    for (int i=mi.MPIlocalCellStart[0]-mi.cellGhostLayerSize; 
         i<mi.MPIlocalCellStart[0]+mi.MPIlocalCellSize[0]+mi.cellGhostLayerSize; i++){

        if (j+size[1]-1 < mi.MPIglobalCellSize[1] // The top of stencil does not exceed maximum Y 
        &&  i+size[0]-1 < mi.MPIglobalCellSize[0] // The right side of stencil does not exceed maximum X
        &&  j > -1                                // The bottom of stencil does not exceed minimum Y
        &&  i > -1                                // The bottom of stencil does not exceed minimum X
        &&  i+size[0]-1 < mi.MPIlocalCellStart[0]+mi.MPIlocalCellSize[0]+mi.cellGhostLayerSize		  
        // The right side of stencil does not exceed mesh partition
        &&  j+size[1]-1 < mi.MPIlocalCellStart[1]+mi.MPIlocalCellSize[1]+mi.cellGhostLayerSize		  
        // The top side of stencil does not exceed mesh partition
        ){
            //if(rank == 1){cout << j << " " << i << endl;}

            indice global();

            interior_.insert(FlatIndic(mi,i,j));
        }

    }}


    return 1;
}
