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
    for (int k=0; k<NS; k++){localbs[k] = -1*localsource[k];}
    for (int k=0; k<Nelem; k++) {localqd[k] = 0.0; localqs[k] = 0.0;}


    PetscFunctionReturn(0);
}

bool Is_Dirichlet(const indice& global){

    return true;
}

// Mark boundary dof in serial
int MarkBndryDOFStokes(bndryVal& bndryStokes, 
                       const MeshInfo& mi, 
                       basis& basis_,
                       BRMixed& br_){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};

        int edge = 0;

        if (Is_Dirichlet(global)){
            // Extract corners coordinates from basis class
            basis_.GetCorners(mi, global);

            vertexSet fullCorners = basis_.corners();

            std::array<int,12> elementDOF = br_.LocalToGlobal(mi,global);

            if (i==0){
                // Count left bottom vertex dof
                // Count left side
                edge = 0; 
            } else if (j==0){
                // Count right bottom vertex dof
                // Count bottom side
                edge = 1;
            } else if (i==mi.MPIglobalCellSize[0]-1){
                // Count right top vertex dof 
                // Count right side
                edge = 2; 
            } else if (j==mi.MPIglobalCellSize[1]-1){
                // Count left top vertex dof
                // Count top side
                edge = 3;
            }

            // Extract two corners representing edge
            vertexSet edgeCorner = {fullCorners.at((edge+3)%4),
                                    fullCorners.at(edge)};
     
            // Extract unit normal vector on the boundary edge
            vertex nu = basis_.unitnormal(edge);

            double len = length(edgeCorner);

            vertex bndryVal = Dirichlet_val(fullCorners.at(edge));

            // x direction
            bndryStokes.insert(std::make_pair<int,bndryInfo>
                               ((int)elementDOF[edge], {edge, bndryVal[0], global}));

            // y direction
            bndryStokes.insert(std::make_pair<int,bndryInfo>
                               ((int)elementDOF[edge+4], {edge+4, bndryVal[1], global}));

            // Assign values to edge supplement bubble function
            double suppVal = AssignBndrySupVal(edgeCorner, nu, gwe, gpe);

            bndryStokes.insert(std::make_pair<int,bndryInfo>
                               ((int)elementDOF[edge+8], {edge+8, suppVal, global}));

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
            // Extract corners of current element to basis functions
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
            } else if (j==mi.MPIglobalCellSize[1]-1){
                // Count top side
                edge = 3;
            }

            // Extract two corners representing edge
            vertexSet edgeCorner = {fullCorners.at((edge+3)%4), 
                                    fullCorners.at(edge)};

            double len = length(edgeCorner);

            // Compute approximated Dirichlet boundary values locally
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
        vertex DiriVal = Dirichlet_val(mapped); 

        l0 += len/2.0*gwe[g]*(DiriVal[0]*vals[0][0] + DiriVal[1]*vals[0][1]);
        l1 += len/2.0*gwe[g]*(DiriVal[0]*vals[1][0] + DiriVal[1]*vals[1][1]);

    }

    work[0] = (a*d-d*d)*(d*l0-b*l1);
    work[1] = (a*d-d*d)*(a*l1-b*l0);

    return work;
}

double AssignBndrySupVal(const vertexSet& edgeCorner,
                         const vertex& nu,
                         const valarray<double>& gwe,
                         const valarray<double>& gpe){

    // Assign value to the degree of freedom of 
    // the supplemental function on the edge
    double work = 0.0;

    // Extract values on both ends of the target edge
    vertex DiriValL = Dirichlet_val(edgeCorner[0]); 
    vertex DiriValR = Dirichlet_val(edgeCorner[1]);

    // Calculate averaged unit normal component
    // of assigned dirichlet boundary values
    double averaged = 0.0;
    for (int g=0; g<gwe.size(); g++) {
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edgeCorner);

        vertex DiriVal = Dirichlet_val(mapped);
        averaged += 1.0/2.0 *gwe[g] *(DiriVal[0] *nu[0] + DiriVal[1]*nu[1]);
    }

    work = averaged - 0.5*((DiriValL[0]+DiriValR[0])*nu[0] + 
                           (DiriValL[1]+DiriValR[1])*nu[1]); 

    work *= 3.0/2.0;

    return work;
}

//! Create Matrix Kg and g for Dirichlet boundary conditions in parallel 
/*
PetscErrorCode CreateDirichletMatVecParallel(Vec * localg,
                                             const bndryVal& bndryvals){

    // Copy vector and matrix
    Vec lg = *localg;   

    // Create with different size
    PetscCall(VecCreate(PETSC_COMM_WORLD, &lg));
    PetscCall(VecSetSizes());

    // Create 
    for (auto & it: bndryvals){


    }


    return PETSC_SUCCESS;
}
*/

//! Create reduced system from full system
PetscErrorCode CreateReducedSystemSerial(ReducedSys * reducedsys,
                                         Mat * fullM,
                                         const bndryVal& bndryval){

    // Copy precalculated full matrix
    Mat fM = *fullM;

    // Get global number of rows and columns from full matrix
    int rows;
    int cols;

    PetscCall(MatGetSize(fM,&rows,&cols));

    assert(rows == cols);

    int bndrySize = (int)bndryval.size();

    int reducedSize = rows - bndrySize;

    // Create reduced system
    PetscCall(MatCreate(PETSC_COMM_WORLD, &(*reducedsys).M));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &(*reducedsys).Kg));

    PetscCall(MatSetSizes(reducedsys->M, PETSC_DECIDE, PETSC_DECIDE, 
              reducedSize, reducedSize);

    PetscCall(MatSetSizes(reducedsys->Kg, PETSC_DECIDE, PETSC_DECIDE, 
              reducedSize, bndrySize);

    PetscCall(MatSetUp(reducedsys->M));
    PetscCall(MatSetUp(reducedsys->Kg));

    // Craete reduced system
	 // Create bndry index array
    std::array<int,bndrySize> bndryIndex;

    int indexg = 0;

    for (auto & it : bndryval){

        bndryIndex[indexg] = (int)it->first;

        bndryInfo info = it->second;

        arrayg[indexg] = it->second

        indexg ++;
    }

    return PETSC_SUCCESS;
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
