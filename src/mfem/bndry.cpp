#include "bndry.h"
#include "myFunc.h"

bool Is_Dirichlet(const indice& global){

    return true;
}

// Mark boundary dof in serial
int MarkBndryDOFStokes(bndryVal& bndryStokes, 
                       const MeshInfo& mi, 
                       basis& basis_,
                       BRMixed& br_){

    // Mark all boundary degree of freedoms
    // Including edge dofs and nodal dofs

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};

        vector<int> edges;

        if (Is_Dirichlet(global)){
            // Extract corners coordinates from basis class
            basis_.GetCorners(mi, global);

            vertexSet fullCorners = basis_.corners();

            std::array<int,12> elementDOF = br_.LocalToGlobal(mi,global);

            if (i==0){
                // Count left bottom vertex dof
                // Count left side
                edges.push_back(0); 
            } 
            if (j==0){
                // Count right bottom vertex dof
                // Count bottom side
                edges.push_back(1);
            } 
            if (i==mi.MPIglobalCellSize[0]-1){
                // Count right top vertex dof 
                // Count right side
                edges.push_back(2); 
            }
            if (j==mi.MPIglobalCellSize[1]-1){
                // Count left top vertex dof
                // Count top side
                edges.push_back(3);
            }

            for (const auto& edge: edges){
                // For each edge 
                // Assign values to only one nodal dofs and one edge dofs
                // Associated local dof are 
                // i, i + 4, i + 8
                // All dofs will be counted without repeating
                vertexSet edgeCorners = {fullCorners.at((edge+3)%4), 
                                         fullCorners.at(edge)};

                vertex nu = basis_.unitnormal(edge);
         
                // Compute values at supplemental bubble function
                double supVal = AssignBndrySupVal(edgeCorners, nu, gwe, gpe);

                // Compute values at nodal dof
                vertex bndryVal = Dirichlet_val(edgeCorners[1]);

                bndryStokes.insert(std::make_pair<int, bndryInfo>
                                   ((int)elementDOF[edge], {edge, bndryVal[0], global}));

                bndryStokes.insert(std::make_pair<int, bndryInfo>
                                   ((int)elementDOF[edge+4], {edge+4, bndryVal[1], global}));

                bndryStokes.insert(std::make_pair<int, bndryInfo>
                                   ((int)elementDOF[edge+8], {edge+8, supVal, global}));
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

        vector<int> edges;

        if (Is_Dirichlet(global)){
            // Extract corners of current element to basis functions
            basis_.GetCorners(mi, global);

            vertexSet fullCorners = basis_.corners();

            std::array<int,8> elementDOF = hdiv_.LocalToGlobal(mi,global);

            if (i==0){
                // Count left side
                edges.push_back(0); 
            }

            if (j==0){
                // Count bottom side
                edges.push_back(1);
            }

            if (i==mi.MPIglobalCellSize[0]-1){
                // Count right side
                edges.push_back(2);
            } 

            if (j==mi.MPIglobalCellSize[1]-1){
                // Count top side
                edges.push_back(3);
            }

            for (const auto& edge : edges){
             //   cout << edge << " " ;

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

            }

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

    vertex nu = basis_.unitnormal(edge);

    for (int g=0; g<gwe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edgeCorner);

        std::array<vertex, 2> vals = hdiv_.ComputeHdivmixed(basis_, mapped, edge); 

        a += len/2.0*gwe[g]*(vals[0][0]*vals[0][0]*nu[0]*nu[0] + 
                             vals[0][1]*vals[0][1]*nu[1]*nu[1]);
        b += len/2.0*gwe[g]*(vals[0][0]*vals[1][0]*nu[0]*nu[0] + 
                             vals[0][1]*vals[1][1]*nu[1]*nu[1]);
        d += len/2.0*gwe[g]*(vals[1][0]*vals[1][0]*nu[0]*nu[0] + 
                             vals[1][1]*vals[1][1]*nu[1]*nu[1]);

        // Get local Dirichlet vector value
        vertex DiriVal = Dirichlet_val(mapped); 

        l0 += len/2.0*gwe[g]*(DiriVal[0]*vals[0][0]*nu[0]*nu[0] + 
                              DiriVal[1]*vals[0][1]*nu[1]*nu[1]);
        l1 += len/2.0*gwe[g]*(DiriVal[0]*vals[1][0]*nu[0]*nu[0] + 
                              DiriVal[1]*vals[1][1]*nu[1]*nu[1]);

    }

    //cout << "Element " << global[0] << " " << global[1] << " l0,l1 : " << l0 << " " << l1 << endl;

    work[0] = (d*l0-b*l1)/(a*d-b*b);
    work[1] = (a*l1-b*l0)/(a*d-b*b);

    return work;
}

double AssignBndrySupVal(const vertexSet& edgeCorner,
                         const vertex& nu,
                         const valarray<double>& gwe,
                         const valarray<double>& gpe){

    // Assign value to the degree of freedom of 
    // the supplemental function on the edge
    // Assign this value to the edge dofs
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
PetscErrorCode CreateReducedSerial(ReducedSys * reducedsys,
                                   Mat * fullM, Mat * fullB,
                                   Vec * fullSource,
                                   const bndryVal& bndryval){

    // Copy precalculated full matrix
    Mat fM = *fullM;
  
    Mat fB = *fullB; 

    Vec fs = *fullSource;

    // Get global number of rows and columns from full matrix
    int rows, cols;

    PetscCall(MatGetSize(fM,&rows,&cols));

    assert(rows == cols);

    int bndrySize = (int)bndryval.size();

    int reducedSize = rows - bndrySize;

    int rowsB, colsB;
    PetscCall(MatGetSize(fB, &rowsB, &colsB));

    assert(colsB < rowsB);
    assert(rowsB == cols);

    // Create reduced system
    PetscCall(MatCreate(PETSC_COMM_WORLD, &reducedsys->M));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &reducedsys->Kg));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &reducedsys->B));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &reducedsys->Bg));

    PetscCall(MatSetSizes(reducedsys->M, PETSC_DECIDE, PETSC_DECIDE, 
              reducedSize, reducedSize));
    PetscCall(MatSetSizes(reducedsys->Kg, PETSC_DECIDE, PETSC_DECIDE, 
              reducedSize, bndrySize));
    PetscCall(MatSetSizes(reducedsys->B, PETSC_DECIDE, PETSC_DECIDE, 
              reducedSize, colsB));
    PetscCall(MatSetSizes(reducedsys->Bg, PETSC_DECIDE, PETSC_DECIDE,
              bndrySize, colsB));

    PetscCall(MatSetUp(reducedsys->M));
    PetscCall(MatSetUp(reducedsys->Kg));
    PetscCall(MatSetUp(reducedsys->B));
    PetscCall(MatSetUp(reducedsys->Bg));

    // Create boundary vector
    PetscCall(VecCreate(PETSC_COMM_WORLD, &reducedsys->g));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &reducedsys->source)); 
    PetscCall(VecSetSizes(reducedsys->g, PETSC_DECIDE, bndrySize));
    PetscCall(VecSetSizes(reducedsys->source, PETSC_DECIDE, reducedSize));
    PetscCall(VecSetUp(reducedsys->g));
    PetscCall(VecSetUp(reducedsys->source));

    double * arrayfullsource;

    PetscCall(VecGetArray(fs, &arrayfullsource));

    int reducedRowIndex = 0;
    int countBndry = 0;
    // Craete reduced system
    for (int row = 0; row < rows; row++){
  
        int bndryIndex = 0;
        int intrIndex  = 0;

        auto itFindRow = bndryval.find(row);
        if (itFindRow == bndryval.end()){
            // boundary vale is empty
            // this dof is interior. Skip all the boundary dofs

            for (int col = 0; col < cols; col++){

                const int idxm = reducedRowIndex;

                // Extract (row, col) value from the pre defined full matrix
                double val;
                MatGetValue(fM, row, col, &val);

                // assign to a const variable.
                const double assignVal = val;

                auto itFind = bndryval.find(col);
                if (itFind != bndryval.end()){
                    // boundary value not empty. 
                    // Means this dof is right on the boundary
                    // 1. Assign boundary value to right hand side
                    // 2. Extract boundary dof related elements from full system

                    const int idxn = bndryIndex; 

                    MatSetValues(reducedsys->Kg, 1, &idxm, 1, &idxn, &assignVal, INSERT_VALUES);

                    bndryIndex ++;                    
                } else {
                    // boundry value is empty
                    // this current dof is interior dof
                    // assign the value to reduced system M

                    const int idxn = intrIndex;

                    MatSetValues(reducedsys->M, 1, &idxm, 1, &idxn, &assignVal, INSERT_VALUES);

                    intrIndex ++;
                }

            }
            double val = arrayfullsource[row];
            VecSetValues(reducedsys->source, 1, &reducedRowIndex, &val , INSERT_VALUES);

            // increment of the reducedRowIndex
            reducedRowIndex ++;
        } else {
            // This dof is on the boundary
            // Put the boundary value into g
            const bndryInfo& tmp = bndryval.at(row);
            VecSetValues(reducedsys->g,1,&countBndry,
                         &tmp.DirichletVal,INSERT_VALUES);
            countBndry ++;
        }
    }

    PetscCall(VecRestoreArray(fs,&arrayfullsource));

    PetscCall(VecAssemblyBegin(reducedsys->g));
    PetscCall(VecAssemblyEnd(reducedsys->g));

    PetscCall(MatAssemblyBegin(reducedsys->M, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(reducedsys->M, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(reducedsys->Kg, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(reducedsys->Kg, MAT_FINAL_ASSEMBLY));

    // Create reduced B
    for (int col =0; col < colsB; col++){
        int bndryIndex = 0;
        int intrIndex = 0;

        for (int row =0; row < rowsB; row++){
            const int idxn = col;
            double val;
            MatGetValue(fB, row, col, &val);
            const double assignVal = val;
 
            auto itFind = bndryval.find(row);
            if (itFind != bndryval.end()){
                const int idxm = bndryIndex;
                MatSetValues(reducedsys->Bg, 1, &idxm, 1, &idxn, &assignVal, INSERT_VALUES);
                bndryIndex ++;
            } else {
                const int idxm = intrIndex;
                MatSetValues(reducedsys->B , 1, &idxm, 1, &idxn, &assignVal, INSERT_VALUES);
                intrIndex ++;
            }
        }

    }

    PetscCall(MatAssemblyBegin(reducedsys->B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(reducedsys->B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(reducedsys->Bg, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(reducedsys->Bg, MAT_FINAL_ASSEMBLY));

    return PETSC_SUCCESS;
}
