#include "passemble.h"

int localDOFPrepare(){
    // Make sure local dof is non overlapping

    return 0;
}

// Assemble sparse matrix parallelly
// Parallel assemble need boundary condition pre allocated
// It is not convenient to delete boundary related row and cols later
PetscErrorCode ParallelMatrixAssembleBlock(const MeshInfo& mi,
                                           basis& basis_,
                                           Hdivmixed& hdiv_,
                                           BRMixed& br_,
                                           PhysProperty * pp, 
                                           System * system){

    PetscFunctionBeginUser;

    // Get gauss points first
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

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
    //PetscCall(MatCreateAij(PETSC_COMM_WORLD, , , br_.getDOF(), br_.getDOF()
    //                                       , , , &(*system).As));

    // Temperatory test
    PetscMPIInt   size; 
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    //PetscCall(MatCreateAij(PETSC_COMM_WORLD, 3, 3, 3*size, 3*size, ));


    return PETSC_SUCCESS;
}
