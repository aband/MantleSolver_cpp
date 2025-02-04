/**!
 * Assemble Darcy-Stokes system
 */

#include "driver.h"

PetscErrorCode Driver::ParallelMatrixAssemble(const Tensor<weights>& allwgtsHD, double ** lHD,
                                              const Tensor<weights>& allwgtsCD, double ** lCD){


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

        // ! Compute cell averaged porosity
        CellAvePorosity(global, allwgtsHD, lHD, allwgtsCD, lCD);

        // ! Compute local values associated to each dofs
        AssignLocMatStokes(global, allwgtsHD, lHD, allwgtsCD, lCD, locmatS);
        AssignLocMatDarcy(global, allwgtsHD, lHD, allwgtsCD, lCD, locmatD);
        AssignLocMatCouple(global, allwgtsHD, lHD, allwgtsCD, lCD, k);

        // ! Load corresponding shape functions
        shape stokesFuncSp(basis_, br_);
        shape darcyFuncSp(basis_, hdiv_);

        // ! Assign local values to global matrix
        if (elemOnBndry(mi, global)){
            // ! Dealing wiht boundary dofs
            AssignLocRedSys(reducedStokes_, locmatS, refArrayStokesEssen_, 
                            mi, bndryStokesEssen_, global, stokesFuncSp, parameter); 
            AssignLocRedSys(reducedDarcy_, locmatD, refArrayDarcyEssen_,
                            mi, bndryDarcyEssen_, global, darcyFuncSp, parameter);
        } else {
            AssignLocRedSys(reducedStokes_, locmatS, refArrayStokesEssen_, mi, global, *br_);
            AssignLocRedSys(reducedDarcy_, locmatD, refArrayDarcyEssen_, mi, global, *hdiv_);
        }

        // Assign coupling K matrix and two C matrices
        // const pressure space not affected by boundary dofs
        PetscCall(MatSetValue(K,nElem,nElem,k,ADD_VALUES));
        PetscCall(MatSetValue(reducedStokes_->C, nElem, nElem, locmatS->C,ADD_VALUES));
        PetscCall(MatSetValue(reducedDarcy_->C, nElem, nElem, locmatD->C, ADD_VALUES));

    }}

    AssembleReducedSys(reducedStokes_);
    AssembleReducedSys(reducedDarcy_);

    PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));

    PetscFunctionReturn(PETSC_SUCCESS);
}
