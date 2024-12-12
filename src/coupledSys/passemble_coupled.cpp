#define COUPLED
#include "passemble.h"

// Assemble sparse matrix parallelly
// Parallel assemble need boundary condition pre allocated
// It is not convenient to delete boundary related row and cols later
// Unlike sequential assembly, we skip assembling full linear system
// We assemble two reduced system directly
PetscErrorCode ParallelMatrixAssemble(const MeshInfo& mi,
                                      basis& basis_,
                                      Phase * phase,
                                      const bndryVal& bndryEssenStokes,
                                      ReducedSys * redsysStokes,
                                      const bndryVal& bndryEssenDarcy,
                                      ReducedSys * redsysDarcy,
                                      Mat * K,
                                      BRMixed& br_,
                                      Hdivmixed& hdiv_,
                                      MLWENO::MLWENOUse * mluse,
                                      int * refArrayStokes, 
                                      int * refArrayDarcy,
                                      const int& bndryDOFStokes,
                                      const int& bndryDOFDarcy,
                                      const std::vector<double>& parameter){

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
        CellAvePorosity(mi, basis_, phase, global, gwf, gpf, mluse);

        // ! Compute local values associated to each dofs
        AssignLocMat(mi, br_  , basis_, locmatS, phase, global, mluse, gwe, gpe, gwf, gpf);
        AssignLocMat(mi, hdiv_, basis_, locmatD, phase, global, mluse, gwe, gpe, gwf, gpf);
        AssignLocMat(mi, br_, hdiv_, basis_, &k, phase, global, mluse, gwf, gpf);

        // ! Load corresponding shape functions
        shape stokesFuncSp(&basis_, &br_);
        shape darcyFuncSp(&basis_, &hdiv_);

        // ! Assign local values to global matrix
        if (elemOnBndry(mi, global)){
            // ! Dealing wiht boundary dofs
            AssignLocRedSys(redsysStokes, locmatS, refArrayStokes, 
                            mi, bndryEssenStokes, global, stokesFuncSp,parameter); 
            AssignLocRedSys(redsysDarcy, locmatD, refArrayDarcy,
                            mi, bndryEssenDarcy, global, darcyFuncSp,parameter);
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
