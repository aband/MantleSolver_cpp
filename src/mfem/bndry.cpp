#include "bndry.h"
#include "myFunc.h"

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

bool Is_Dirichlet(const indice& global){

    return true;
}

// Mark boundary dof in serial
int MarkBndryDOFStokes(bndryVal& bndryStokes, 
                       const MeshInfo& mi, BRMixed& br_){

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};
        if (Is_Dirichlet(global)){

            std::array<int,12> elementDOF = br_.LocalToGlobal(mi,global);

            if (i==0){
                // Count left bottom vertex dof
                // Count left sid
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
        } // else (for Neumann situation) 
    }}

    return 0;
}

int MarkBndryDOFDarcy(bndryVal& bndryDarcy, 
                      const MeshInfo& mi, 
                      basis& basis_,
                      Hdivmixed& hdiv_){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};

        int edge = 0;

        if (Is_Dirichlet(global)){
            // Assign corners of current element to basis functions
            basis_.GetCorners(mi, global);

            vertexSet fullCorners = basis_.corners();

            std::array<int,8> elementDOF = hdiv_.LocalToGlobal(mi,global);

            if (i==0){
                // Count left side
                edge = 0; 
            } else if (j==0){
                // Count bottom side
                edge = 1;
            } else if (i==mi.MPIglobalCellSize[0]-1){
                // Count right side
                edge = 2;
            } else if (j=mi.MPIglobalCellSize[0]-1){
                // Count top side
                edge = 3;
            }

            // Extract two corners representing edge
            vertexSet edgeCorner = {fullCorners.at((edge+3)%4), 
                                    fullCorners.at(edge)};

            double len = length(edgeCorner);

            std::array<double, 2> dVals = AssignBndryValsDarcy(global, edge, basis_,hdiv_,
                                                               edgeCorner, len, gwe, gpe);

            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[edge], {edge,dVals[0],global}));

            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                               ((int)elementDOF[edge+4], {edge+4,dVals[1],global}));

        } // else (save later for neumann boundary condition)
    }}

    return 0;
}

int AssignBndryValsStokes(bndryVal& bndryStokes, BRMixed& br_){

    // Assign point wise value directly

    for(auto& it: bndryStokes){

    }

    return 0;
}

std::array<double,2> AssignBndryValsDarcy(const indice& global,
                                          const int& edge,
                                          basis& basis_,
                                          Hdivmixed& hdiv_,
                                          const vertexSet& edgeCorner,
                                          const double& len,
                                          const valarray<double>& gwe,
                                          const valarray<double>& gpe){

    std::array<double, 2> work;

    // Initialize the local linear system variables
    double a = 0, b = 0, d = 0;

    // Initialize the right hand side vector components
    double l0 = 0, l1 = 0;

    // Assign Dirichlet boundary values to Darcy problem
    // requires a L2 projection.
    // Vector based basis function cannot assign Dirichlet 
    // boundary condition directly.
    // A first order approximation minimization L2 error.

    for (int g=0; g<gwe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edgeCorner);
        std::array<vertex, 2> vals = hdiv_.ComputeHdivmixed(basis_, mapped, edge); 
        a += len/2.0*gwe[g]*(vals[0][0]*vals[0][0] + vals[0][1]*vals[0][1]);
        b += len/2.0*gwe[g]*(vals[0][0]*vals[1][0] + vals[0][1]*vals[1][1]);
        d += len/2.0*gwe[g]*(vals[1][0]*vals[1][0] + vals[1][1]*vals[1][1]);

        // Get local Dirichlet vector value
        vertex dVal = Dirichlet_val(mapped); 

        l0 += len/2.0*gwe[g]*(dVal[0]*vals[0][0] + dVal[1]*vals[0][1]);
        l1 += len/2.0*gwe[g]*(dVal[0]*vals[1][0] + dVal[1]*vals[1][1]);

    }

    work[0] = (a*d-d*d)*(d*l0-b*l1);
    work[1] = (a*d-d*d)*(a*l1-b*l0);

    return work;
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
