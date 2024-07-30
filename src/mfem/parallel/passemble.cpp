#include "passemble.h"

PetscErrorCode ParallelAssembleTest(){

    // Temperatory test
    PetscMPIInt   size, rank; 
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    Mat test;
    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, 3, 3, 3*size, 3*size, 2, NULL , 2, NULL, &test));

    PetscCall(MatSetUp(test));

    PetscCall(MatSetValue(test, 0, 0, 1.0, ADD_VALUES));

    int ownerm, ownern;
    PetscCall(MatGetOwnershipRange(test, &ownerm, &ownern));

    PetscPrintf(PETSC_COMM_SELF, "rank = %d, ownership m = %d, ownership n = %d \n",
                rank, ownerm, ownern);

    PetscCall(MatAssemblyBegin(test, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(test, MAT_FINAL_ASSEMBLY));

    return PETSC_SUCCESS;
}

/*
template <typename T>
int CreateRefMap(T& funcSp, int * refArray, 
                 const MeshInfo& mi, int * bndryDOFEssen){

    // Each processor has to create its own mapping
    // Control Essential dof only

    // !!!!!! Caution !!!!!!!
    // This function has not been finished
    // Cannot do natural boundary condition yet

    int bndryIndex = 0;
    int intrIndex  = 0;

    for (int dof=0; dof<funcSp.getDOF(); dof++){

        if (funcSp.onBndry(mi, dof)){

            refArray[dof] = intrIndex;
            intrIndex ++;
        }else {
            refArray[dof] = bndryIndex;
            bndryIndex ++;
        }

    }

    *bndryDOFEssen = intrIndex;

    return 0;  
}
*/

// Assemble sparse matrix parallelly
// Parallel assemble need boundary condition pre allocated
// It is not convenient to delete boundary related row and cols later
// Unlike sequential assembly, we skip assembling full linear system
// We assemble two reduced system directly
PetscErrorCode ParallelMatrixAssemble(const MeshInfo& mi,
                                      basis& basis_,
                                      PhysProperty * pp,
                                      const bndryVal& bndryEssenStokes,
                                      ReducedSys * redsysStokes,
                                      const bndryVal& bndryEssenDarcy,
                                      ReducedSys * redsysDarcy,
                                      Mat * K,
                                      BRMixed& br_,
                                      Hdivmixed& hdiv_,
                                      int * refArrayStokes, 
                                      int * refArrayDarcy,
                                      const int& bndryDOFStokes,
                                      const int& bndryDOFDarcy){

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

//    int bndryDOFStokes = 0.0;
//    int bndryDOFDarcy = 0.0;

//    int * refArrayStokes = new int[br_.getDOF()];
//    int * refArrayDarcy  = new int[hdiv_.getDOF()];

//    CreateRefMap(br_, refArrayStokes, mi, &bndryDOFStokes);
//    CreateRefMap(hdiv_, refArrayDarcy, mi, &bndryDOFDarcy);

    int reducedDOFStokes = br_.getDOF() - bndryDOFStokes;

    int reducedDOFDarcy = hdiv_.getDOF() - bndryDOFDarcy;

    // Create parallel sparse matrix ===============================================
    // MatCreateAIJ(MPI_Comm comm, PetscInt m, PetscInt n, 
    //                             PetscInt M, PetscInt N, 
    //                             PetscInt d_nz, const PetscInt d_nnz[], 
    //                             PetscInt o_nz, const PetscInt o_nnz[], Mat *A)
    // m number of local rows , n = m for square matrix
    // M number of gloabl rows, N number of gloabl columns
    // d_nz number of nonzeros per row in Diagonal portion of local submatrix
    // d_nnz array containing the number of nonzero blocks
    // o_nz, number of nonzero blocks per block
    // o_nnz, array containing the number of nonzero blocks
    // =============================================================================

    PrepareReducedSys(redsysStokes, reducedDOFStokes, bndryDOFStokes, 
                      totalElem, 30, 30, 4, 4);
    PrepareReducedSys(redsysDarcy, reducedDOFDarcy, bndryDOFDarcy, 
                      totalElem, 14, 14, 2, 2);

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           totalElem, totalElem, 
                           1, NULL, 0, NULL, K));
    PetscCall(MatSetUp(*K));

    // =============================================================================

//    int vsize;
//    VecGetSize(redsysStokes->g, &vsize);
//    cout << vsize << endl;
//    cout << bndryDOFStokes << endl;

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
        basis_.GetCorners(mi, global);

        // ! Compute cell averaged porosity
        CellAvePorosity(mi, pp, basis_, gwf, gpf);

        // ! Compute local values associated to each dofs
        AssignLocMat(mi, br_  , basis_, locmatS, pp, gwe, gpe, gwf, gpf);
        AssignLocMat(mi, hdiv_, basis_, locmatD, pp, gwe, gpe, gwf, gpf);
        AssignLocMat(mi, br_, hdiv_, basis_, pp, &k, gwf, gpf);

        // ! Load corresponding shape functions
        shape stokesFuncSp(&basis_, &br_);
        shape darcyFuncSp(&basis_, &hdiv_);

        // ! Assign local values to global matrix
        if (elemOnBndry(mi, global)){
            AssignLocRedSys(redsysStokes, locmatS, refArrayStokes, 
                            mi, bndryEssenStokes, global, stokesFuncSp); 
            AssignLocRedSys(redsysDarcy, locmatD, refArrayDarcy,
                            mi, bndryEssenDarcy, global, darcyFuncSp);
        } else {
            AssignLocRedSys(redsysStokes, locmatS, refArrayStokes, mi, global, br_);
            AssignLocRedSys(redsysDarcy, locmatD, refArrayDarcy, mi, global, hdiv_);
        }

        // Assign coupling K matrix and two C matrices
        // const pressure space not affected by boundary dofs
        PetscCall(MatSetValue(*K,nElem,nElem,k,ADD_VALUES));
        PetscCall(MatSetValue(redsysStokes->C, nElem, nElem, locmatS->C,ADD_VALUES));
        PetscCall(MatSetValue(redsysDarcy->C, nElem, nElem, locmatD->C, ADD_VALUES));

    }}

    AssembleReducedSys(redsysStokes);
    AssembleReducedSys(redsysDarcy);

    PetscCall(MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY));

    return PETSC_SUCCESS;
}
