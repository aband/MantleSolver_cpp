#include "bndry.h"

PetscErrorCode AssignValuesRHS(int NS, int ND, int Nelem,
                               Vec * A, Vec * B,
                               RHSVector * rhsv,
                               const bndryVal& bndryStokes,
                               const bndryVal& bndryDarcy){

    PetscFunctionBeginUser;

    Vec a, b;

    a = *A;
    b = *B;

    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, ND, &rhsv->ad));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, NS, &rhsv->bs));

    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, Nelem, &rhsv->qd));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, Nelem, &rhsv->qs));

    // Assign local arrays
    double *localad;
    double *localbs;
    double *localqd;
    double *localqs;

    double *localsource;

    PetscCall(VecGetArray(rhsv->ad, &localad));
    PetscCall(VecGetArray(rhsv->bs, &localbs));
    PetscCall(VecGetArray(rhsv->qd, &localqd));
    PetscCall(VecGetArray(rhsv->qs, &localqs));

    PetscCall(VecGetArray(rhsv->source, &localsource));

    // Initialize local arrays
    for (int k=0; k<ND; k++){localad[k] = 0.0;}
    for (int k=0; k<NS; k++){localbs[k] = 0.0;}
    for (int k=0; k<Nelem; k++) {localqd[k] = 0.0; localqs[k] = 0.0;}




    PetscFunctionReturn(0);
}

// Mark boundary dof in serial
int MarkBndryDOFStokes(bndryVal& bndryStokes, const MeshInfo& mi, BRMixed& br_){

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};

        std::array<int,12> elementDOF = br_.LocalToGlobal(mi,global);

        if (i==0){
            // Count left bottom vertex dof
            // Count left side
            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[0], {0,0.0,global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[4], {4,0.0,global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[8], {8,0.0,global}));
        } else if (j==0){
            // Count right bottom vertex dof
            // Count bottom side
            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[1], {1,0.0,global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[5], {5,0.0,global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[9], {9,0.0,global}));
 
        } else if (i==mi.MPIglobalCellSize[0]-1){
            // Count right top vertex dof 
            // Count right side
            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[2], {2,0.0,global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[6], {6,0.0,global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[10], {10,0.0,global}));
 
        } else if (j=mi.MPIglobalCellSize[0]-1){
            // Count left top vertex dof
            // Count top side
            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[3], {3,0.0,global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[7], {7,0.0,global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[11], {11,0.0,global}));
 
        }

    }}

    return 0;
}

int MarkBndryDOFDarcy(bndryVal& bndryDarcy, const MeshInfo& mi, Hdivmixed& hdiv_){

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};

        std::array<int,8> elementDOF = hdiv_.LocalToGlobal(mi,global);

        if (i==0){
            // Count left side
            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[0], {0,0.0,global}));

            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[4], {4,0.0,global}));

        } else if (j==0){
            // Count bottom side
            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[1], {1,0.0,global}));

            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[5], {5,0.0,global}));

        } else if (i==mi.MPIglobalCellSize[0]-1){
            // Count right side
            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[2], {2,0.0,global}));

            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[6], {6,0.0,global}));

        } else if (j=mi.MPIglobalCellSize[0]-1){
            // Count top side
            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[3], {3,0.0,global}));

            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[7], {7,0.0,global}));

        }
    }}

    return 0;
}

int ComputeBndryValsStokes(bndryVal& bndryStokes, BRMixed& br_,
                           const valarray<double>& gwe,
                           const valarray<double>& gpe){

    for(auto& it: bndryStokes){

    }

    return 0;
}

int ComputeBndryValsDarcy(bndryVal& bndryDarcy, Hdivmixed& hdiv_,
                          const valarray<double>& gwe,
                          const valarray<double>& gpe){


    return 0;
}

PetscErrorCode CreateRHS(const MeshInfo& mi,
                         basis& basis_,
                         Hdivmixed& hdiv_,
                         BRMixed& br_,
                         PhysProperty * physproperty,
                         RHSVector * rhsv){

    // Create right hand side vector
    PetscFunctionBeginUser;

    // Create auxilliary vectors
//    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, ));


    PetscFunctionReturn(0);
}
