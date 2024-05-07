#include "pbndry.h"

// Mark boundary dofs in parallel computing
int ParallelMarkBndryDOFs(const MeshInfo& mi,
                          bndryVal& bndryDiri,
                          bndryVal& bndryNeum,
                          basis& basis_,
                          PhysProperty * pp,
                          BRMixed& br_){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    int jstart = mi.MPIlocalCellStart[1];
    int istart = mi.MPIlocalCellStart[0];

    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    PetscPrintf(PETSC_COMM_SELF, "Current rank is %d \n", rank);

    // Loop through local partition of physical domain
    // No ghost layer included
    // No overlapping element index across 
    for (int j=jstart; j<jstart+mi.MPIlocalCellSize[1]; j++){
    for (int i=istart; i<istart+mi.MPIlocalCellSize[0]; i++){

        PetscPrintf(PETSC_COMM_SELF, "j = %d, i=%d ",j,i);

    }PetscPrintf(PETSC_COMM_SELF, "\n");}

    return 0;
}
