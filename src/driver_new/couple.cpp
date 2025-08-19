#include "driver.h"

// Functions regarding coupled computation
int Driver::SolveFlow(int maxIter, double tolUzawa, double ** lHD, double ** lCD){

    ParallelMatrixAssemble(lHD, lCD);

    int nelem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    CreateLinearSys(reducedStokes_, nelem);
    CreateLinearSys(reducedDarcy_, nelem);

    CreateCoupledSystem(reducedStokes_, reducedDarcy_, Result_, &K);

    CoupledUzawa(Result_, tolUzawa, maxIter);

    return 1;
}

int Driver::ParallelMatrixAssemble(double ** lHD, double ** lCD){

    PetscMPIInt   size, rank; 
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscFunctionBeginUser;

    // Get gauss points first
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    // Calculate dofs 
    int totalElem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    int reducedDOFStokes = br_->getDOF() - bndryDOFStokes_;

    int reducedDOFDarcy = hdiv_->getDOF() - bndryDOFDarcy_;

    PrepareReducedSys(reducedStokes_, reducedDOFStokes, bndryDOFStokes_, 
                      totalElem, 30, 30, 4, 4);
    PrepareReducedSys(reducedDarcy_, reducedDOFDarcy, bndryDOFDarcy_, 
                      totalElem, 14, 14, 2, 2);

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           totalElem, totalElem, 
                           1, NULL, 0, NULL, &K));
    PetscCall(MatSetUp(K));

    // ===================================================================

    LocMat * locmatS = new LocMat;
    LocMat * locmatD = new LocMat;

    double k = 0.0;

    // ! Loop local portion of physical domain
    int istart = mi.MPIlocalCellStart[0];
    int jstart = mi.MPIlocalCellStart[1];

    for (int j=jstart; j<jstart + mi.MPIlocalCellSize[1]; j++){
    for (int i=istart; i<istart + mi.MPIlocalCellSize[0]; i++){

        // ! Get global element index
        indice global {i,j};

        int nElem = FlatIndic(mi,global);

        // ! Extract corners of this element
        basis_->GetCorners(mi, global);




    }}

    AssembleReducedSys(reducedStokes_);
    AssembleReducedSys(reducedDarcy_);

    PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));

    PetscFunctionReturn(PETSC_SUCCESS);

    return 1;
}
