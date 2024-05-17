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

    PetscCall(MatView(test, PETSC_VIEWER_STDOUT_WORLD));

    return PETSC_SUCCESS;
}

inline bool elemOnBndry(const MeshInfo& mi,
                        const indice& global){
     
    if (global[0] == 0 ||
        global[1] == 0 ||
        global[0] == mi.MPIglobalCellSize[0] - 1 ||
        global[1] == mi.MPIglobalCellSize[1] - 1){

        return true;

    } else {

        return false;

    }
}

template <typename T>
inline int CreateRefMap(T& funcSp, int * refArray, 
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

inline int PrepareReducedSys(ReducedSys * redsys, 
                             int reducedDOF, int bndrySize, int totalElem,
                             int Adnz, int Aonz,
                             int Bdnz, int Bonz){

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             reducedDOF, reducedDOF, 
                                             Adnz, NULL, Aonz, NULL, &redsys->M));
    PetscCall(MatSetUp(redsys->M));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             reducedDOF, bndrySize,
                                             Adnz, NULL, Aonz, NULL, &redsys->Kg));  
    PetscCall(MatSetUp(redsys->Kg));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             reducedDOF, totalElem,
                                             Bdnz, NULL, Bonz, NULL, &redsys->B));   
    PetscCall(MatSetUp(redsys->B));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             bndrySize, totalElem,
                                             Bdnz, NULL, Bonz, NULL, &redsys->Bg));  
    PetscCall(MatSetUp(redsys->Bg));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             totalElem, totalElem,
                                             1, NULL, 0, NULL, &redsys->C));   
    PetscCall(MatSetUp(redsys->C));

    // Create Corresponding vectors
    // Get ownership first
    // g vector has the size of bndrySize
    // source vector has the size of reducedDOF 
    int m, n; 
    PetscCall(MatGetOwnershipRange(redsys->Bg, &m, &n));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys->g));

    PetscCall(MatGetOwnershipRange(redsys->B, &m, &n));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys->source));

    return 0;
}

inline int AssembleReducedSys(ReducedSys * redsys){

    PetscCall(MatAssemblyBegin(redsys->M, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->M, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(redsys->Kg, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->Kg, MAT_FINAL_ASSEMBLY));

    PetscCall(MatAssemblyBegin(redsys->B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(redsys->Bg, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->Bg, MAT_FINAL_ASSEMBLY));

    PetscCall(MatAssemblyBegin(redsys->C, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->C, MAT_FINAL_ASSEMBLY));

    PetscCall(VecAssemblyBegin(redsys->g));
    PetscCall(VecAssemblyEnd(redsys->g));

    PetscCall(VecAssemblyBegin(redsys->source));
    PetscCall(VecAssemblyEnd(redsys->source));

    return 0;
}

// Used for elements on the boundary (has dofs on the boundary)
template <typename T> 
inline int AssignLocRedSys(ReducedSys * redsys,
                           LocMat * loc,
                           int * ref,
                           const MeshInfo& mi,
                           const bndryVal& bndryEssen,
                           const indice& global,
                           shape<T>& funcSp){

    // ! Get global index of local dofs 
    const std::vector<int> elemDofs = funcSp.LocalToGlobal(mi, global);

    const int idxn = FlatIndic(mi, global);

    for (int row=0; row<elemDofs.size(); row++){
        // View it as the row index
        const int idxm = ref[elemDofs.at(row)];
        const double valB = loc->B.at(row);

        if (funcSp.onBndry(mi, elemDofs.at(row))){
            // This dof is on the boundary
            // Should be assign to Bg
            PetscCall(MatSetValues(redsys->Bg, 1, &idxm, 1, &idxn, &valB, 
                                   ADD_VALUES));

            // At the same time insert essential boundary value to rhs vector
            auto itFind = bndryEssen.find(elemDofs.at(row));
            if (itFind != bndryEssen.end()){
                const bndryInfo& tmp = bndryEssen.at(elemDofs.at(row));
                PetscCall(VecSetValues(redsys->g, 1, &idxm, &tmp.val, INSERT_VALUES));
            }
        } else {
            // This dof is not on the boundary
            // Should be assigned to B instead
            PetscCall(MatSetValues(redsys->B , 1, &idxm, 1, &idxn, &valB, 
                                   ADD_VALUES));

            // This dof is not on the boundary
            // This dof will contribute to source term
            double vals = loc->f.at(row);
            PetscCall(VecSetValues(redsys->source, 1, &idxm, &vals, ADD_VALUES));

            for (int col=0; col<elemDofs.size(); col++){
                const int cidxn  = ref[elemDofs.at(col)];
                const double val = loc->A.at(row+col*elemDofs.size()); 

                if (funcSp.onBndry(mi,elemDofs.at(col))){
                    // It is a non bndry dof - bndry dof interaction
                    // Val assigned to M
                    PetscCall(MatSetValues(redsys->Kg, 1, &idxm, 1, &cidxn, 
                                           &val, ADD_VALUES));
                } else {
                    // It is a non bndry dof - non bndry dof interaction
                    // Val assigned to M
                    PetscCall(MatSetValues(redsys->M, 1, &idxm, 1, &cidxn, 
                                           &val, ADD_VALUES));
                }
            }
        }
    }

    return 0;
}

// Used for interior elements (no need to identify boundary dofs)
// Set multiple values at the same time
// need const int id array
// Stokes and Darcy part are separated
inline int AssignLocRedSys(ReducedSys * redsys,
                           LocMat * loc,
                           int * ref,
                           const MeshInfo& mi,
                           const indice& global,
                           Hdivmixed& hdiv_){

    const int idxm = FlatIndic(mi, global);

    std::array<int, 8> locDof = hdiv_.LocalToGlobal(mi, global);

    const int IS[8] = {ref[locDof[0]], ref[locDof[1]], 
                       ref[locDof[2]], ref[locDof[3]],
                       ref[locDof[4]], ref[locDof[5]], 
                       ref[locDof[6]], ref[locDof[7]]};

    const double locB[8] = {loc->B.at(0), loc->B.at(1),
                            loc->B.at(2), loc->B.at(3),
                            loc->B.at(4), loc->B.at(5),
                            loc->B.at(6), loc->B.at(7)};

    // Only fill B and M matrix
    PetscCall(MatSetValues(redsys->B, 8, IS, 1, &idxm, locB, ADD_VALUES));

    for (unsigned int l=0; l<8; l++){
        const double locA[8] = {loc->A.at(0+l*8), loc->A.at(1+l*8),
                                loc->A.at(2+l*8), loc->A.at(3+l*8),
                                loc->A.at(4+l*8), loc->A.at(5+l*8),
                                loc->A.at(6+l*8), loc->A.at(7+l*8)};
        const int Aidxm = IS[l];
        PetscCall(MatSetValues(redsys->M,1,&Aidxm,8,IS,locA,ADD_VALUES));
    }

    // Fill Source vector
    const double locf[8] = {loc->f.at(0), loc->f.at(1),
                            loc->f.at(2), loc->f.at(3),
                            loc->f.at(4), loc->f.at(5),
                            loc->f.at(6), loc->f.at(7)};

    PetscCall(VecSetValues(redsys->source, 8, IS, locf, ADD_VALUES));

    return 0;
}

inline int AssignLocRedSys(ReducedSys * redsys,
                           LocMat * loc,
                           int * ref,
                           const MeshInfo& mi,
                           const indice& global,
                           BRMixed& br_){

    const int idxm = FlatIndic(mi, global);

    std::array<int, 12> locDof = br_.LocalToGlobal(mi, global);

    const int IS[12] = {ref[locDof[0]], ref[locDof[1]], 
                        ref[locDof[2]], ref[locDof[3]],
                        ref[locDof[4]], ref[locDof[5]], 
                        ref[locDof[6]], ref[locDof[7]],
                        ref[locDof[8]], ref[locDof[9]],
                        ref[locDof[10]], ref[locDof[11]]};

    const double locB[12] = {loc->B.at(0), loc->B.at(1),
                             loc->B.at(2), loc->B.at(3),
                             loc->B.at(4), loc->B.at(5),
                             loc->B.at(6), loc->B.at(7),
                             loc->B.at(8), loc->B.at(9),
                             loc->B.at(10), loc->B.at(11)};

    // Only fill B and M matrix
    PetscCall(MatSetValues(redsys->B, 12, IS, 1, &idxm, locB, ADD_VALUES));

    for (unsigned int l=0; l<12; l++){
        const double locA[12] = {loc->A.at(0+l*12), loc->A.at(1+l*12),
                                 loc->A.at(2+l*12), loc->A.at(3+l*12),
                                 loc->A.at(4+l*12), loc->A.at(5+l*12),
                                 loc->A.at(6+l*12), loc->A.at(7+l*12),
                                 loc->A.at(8+l*12), loc->A.at(9+l*12),
                                 loc->A.at(10+l*12), loc->A.at(11+l*12)};
        const int Aidxm = IS[l];
        PetscCall(MatSetValues(redsys->M,1,&Aidxm,12,IS,locA,ADD_VALUES));
    }

    // Only contribute to source vector not g vector

    const double locf[12] = {loc->f.at(0), loc->f.at(1),
                             loc->f.at(2), loc->f.at(3),
                             loc->f.at(4), loc->f.at(5),
                             loc->f.at(6), loc->f.at(7),
                             loc->f.at(8), loc->f.at(9),
                             loc->f.at(10), loc->f.at(11)};

    PetscCall(VecSetValues(redsys->source, 12, IS, locf, ADD_VALUES));

    return 0;
}

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
                                      Hdivmixed& hdiv_){

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

    int bndryDOFStokes = 0.0;
    int bndryDOFDarcy = 0.0;

    int * refArrayStokes = new int[br_.getDOF()];
    int * refArrayDarcy  = new int[hdiv_.getDOF()];

    CreateRefMap(br_, refArrayStokes, mi, &bndryDOFStokes);
    CreateRefMap(hdiv_, refArrayDarcy, mi, &bndryDOFDarcy);

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
                      totalElem, 30, 22, 4, 4);
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

    free(refArrayStokes);
    free(refArrayDarcy);

    return PETSC_SUCCESS;
}
