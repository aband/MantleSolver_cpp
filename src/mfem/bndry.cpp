#include "bndry.h"
#include "myFunc.h"

void markBndryEdge(const MeshInfo& mi,
                   vector<int>& edges,
                   const int& i, const int& j){

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
}

/*
bool isCorner(const MeshInfo& mi, 
              const indice& global){
    // Judge if the element is a corner element or not 
    if (global[0] == 0 || global[0] == mi.MPIglobalCellSize[0]-1){
        if (global[1] == 0 || global[1] = mi.MPIglobalCellSize[1]-1){
            return true;
        }
    } else {
        return false;
    }
}
*/

// Mark boundary dof in a general way
// Dirichlet and Neumann boundary condition
int MarkBndryDOFStokes(bndryVal& bndryDiri,
                       bndryVal& bndryNeum,
                       const MeshInfo& mi,
                       basis& basis_,
                       BRMixed& br_,
                       PhysProperty * pp){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        // Global element index
        indice global{i,j};

        // Extract corners of this element
        basis_.GetCorners(mi, global);

        vertexSet fullCorners = basis_.corners();

        // Get global numbering of the dofs associating with this element
        std::array<int, 12> elementDOF = br_.LocalToGlobal(mi, global);

        // Mark all the edges of this element that laying on the boundary
        vector<int> edges;

        markBndryEdge(mi, edges, i, j); 

        for (const auto& edge: edges){
            // Get corners corresponding to this boundary edge

            vertexSet edgeCorners = {fullCorners.at((edge+3)%4),
                                     fullCorners.at(edge)};

            // Get unit normal vector to this boundary edge
            vertex nu = basis_.unitnormal(edge);

            // Get dirichlet boundary nodal value
            // and supplemental bubble function value
            vertex bndryVal = bndryVs(edgeCorners[1], pp);
            double supVal = AssignBndrySupVal(edgeCorners, nu, gwe, gpe, pp);

            std::array<double,3> tmpVal {bndryVal[0], bndryVal[1], supVal};

            double neumVal = 0.0;

            // Three dofs associated with this edge are counted here
            for (int dofi = 0; dofi < 3; dofi++){
                int locdof = edge + dofi*4;

                switch (bndryTypeMarker(mi, global, locdof)){
                    case dirichlet:
                        // Dirichlet boundary condition
                        bndryDiri.insert(std::make_pair<int, bndryInfo>
                             ((int)elementDOF[locdof],{locdof, tmpVal[dofi], global}));

                        break;
                     
                    case neumann:
                        // Neumann boundary condition
                        // Calculate boundry integration relates to this dof
                        neumVal = neumValStokes(
                                  mi,global,edge,locdof,dofi,basis_,br_,gwe,gpe,pp);
                        bndryNeum.insert(std::make_pair<int, bndryInfo>
                             ((int)elementDOF[locdof],{locdof, neumVal, global}));
                        break;

                    case missed:
                        cout << "This dof is missed." << endl;
                        break;
                }
            }
        }
    }}

    return 0; 
}

// Mark boundary dof in serial
// Dirichlet only function
int MarkBndryDOFStokes(bndryVal& bndryStokes, 
                       const MeshInfo& mi, 
                       basis& basis_,
                       BRMixed& br_,
                       PhysProperty * pp){

    // Mark all boundary degree of freedoms
    // Including edge dofs and nodal dofs

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};

        vector<int> edges;

        // Extract corners coordinates from basis class
        basis_.GetCorners(mi, global);

        vertexSet fullCorners = basis_.corners();

        std::array<int,12> elementDOF = br_.LocalToGlobal(mi,global);

        markBndryEdge(mi, edges, i, j);

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
            double supVal = AssignBndrySupVal(edgeCorners, nu, gwe, gpe, pp);

            // Compute values at nodal dof
            //vertex bndryVal = Dirichlet_val(edgeCorners[1]);
            vertex bndryVal = bndryVs(edgeCorners[1], pp); 

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                 ((int)elementDOF[edge], {edge, bndryVal[0], global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                 ((int)elementDOF[edge+4], {edge+4, bndryVal[1], global}));

            bndryStokes.insert(std::make_pair<int, bndryInfo>
                 ((int)elementDOF[edge+8], {edge+8, supVal, global}));
        }
    }}

    return 0;
}

// Dirichlet only function
int MarkBndryDOFDarcy(bndryVal& bndryDarcy, 
                      const MeshInfo& mi, 
                      basis& basis_,
                      Hdivmixed& hdiv_,
                      PhysProperty * pp){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};

        vector<int> edges;

        // Extract corners of current element to basis functions
        basis_.GetCorners(mi, global);

        vertexSet fullCorners = basis_.corners();

        std::array<int,8> elementDOF = hdiv_.LocalToGlobal(mi,global);

        markBndryEdge(mi, edges, i, j);

        for (const auto& edge : edges){

            // Extract two corners representing edge
            vertexSet edgeCorner = {fullCorners.at((edge+3)%4), 
                                    fullCorners.at(edge)};

            double len = length(edgeCorner);

            // Compute approximated Dirichlet boundary values locally
            std::array<double, 2> dVals = AssignBndryValsDarcy(global, edge, 
                 basis_,hdiv_, pp, edgeCorner, len, gwe, gpe);
    
            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                    ((int)elementDOF[edge], {edge,dVals[0],global}));

            bndryDarcy.insert(std::make_pair<int, bndryInfo>
                    ((int)elementDOF[edge+4], {edge+4,dVals[1],global}));
        }
    }}

    return 0;
}

std::array<double,2> AssignBndryValsDarcy(const indice& global,
                                          const int& edge,
                                          basis& basis_,
                                          Hdivmixed& hdiv_,
                                          PhysProperty * pp, 
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
        //vertex DiriVal = Dirichlet_val(mapped); 
        vertex DiriVal = bndryu(mapped,pp);

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
                         const valarray<double>& gpe,
                         PhysProperty * pp){

    // Assign value to the degree of freedom of 
    // the supplemental function on the edge
    // Assign this value to the edge dofs
    double work = 0.0;

    // Extract values on both ends of the target edge
    //vertex DiriValL = Dirichlet_val(edgeCorner[0]); 
    //vertex DiriValR = Dirichlet_val(edgeCorner[1]);
    vertex DiriValL = bndryVs(edgeCorner[0], pp); 
    vertex DiriValR = bndryVs(edgeCorner[1], pp);

    // Calculate averaged unit normal component
    // of assigned dirichlet boundary values
    double averaged = 0.0;
    for (int g=0; g<gwe.size(); g++) {
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edgeCorner);

        //vertex DiriVal = Dirichlet_val(mapped);
        vertex DiriVal = bndryVs(mapped, pp);
        averaged += 1.0/2.0 *gwe[g] *(DiriVal[0] *nu[0] + DiriVal[1]*nu[1]);
    }

    work = averaged - 0.5*((DiriValL[0]+DiriValR[0])*nu[0] + 
                           (DiriValL[1]+DiriValR[1])*nu[1]); 

    work *= 3.0/2.0;

    if ((nu[0]+nu[1])<0){
        work *= -1;
    }

    return work;
}

double neumValStokes(const MeshInfo& mi,
                     const indice& global, 
                     const int& edge,
                     const int& local,
                     const int& dofi, 
                     basis&   basis_,
                     BRMixed& br_,
                     const valarray<double>& gwe,
                     const valarray<double>& gpe,
                     PhysProperty * pp){

    double work = 0.0;

    // Owner element of the shape function
    vector<indice> owner;
    vector<int> localdof;
    vector<int> intEdges;

    // Add the first integral range to the vector
    owner.push_back(global);
    localdof.push_back(local);
    intEdges.push_back(edge);
 
    bool is_corner = false;

    indiceSet addSet {{0,-1},{1,0},{0,1},{-1,0}};

    // Find the second integral range
    if (local < 8){
        // Nodal dof

        switch (edge){
            case 0:
                // This element is on the left boundary
                // need edge 0 (i,j) and (i,j-1), expect corner
                if (global[1] == 0 ){ // corner
                    is_corner = true; 
                } else { // edge
                    is_corner = false;
               }

                break;
            case 1:
                // This element is on the bottom boundary
                // need edge 1 (i,j) and (i+1,j), expect corner
                if (global[0] == mi.MPIglobalCellSize[0]-1){ // corner
                    is_corner = true;
                } else {
                    is_corner = false;
                }

                break;
            case 2:
                // This element is on the right boundary
                // need edge 2 (i,j) and (i,j+1), expect corner
                if (global[1] == mi.MPIglobalCellSize[1]-1){
                    is_corner = true;
                } else {
                    is_corner = false;
                }

                break;
            case 3:
                // This element is on the top boundary
                // need edge 3 (i,j) and (i-1,j), expect corner
                if (global[1] == 0){
                    is_corner = true;
                } else {
                    is_corner = false;
                }

                break;
            default:
                cout << "Undefined edge." << endl;
                break;
        }

        if (is_corner){
            owner.push_back(global);
            intEdges.push_back((edge+1)%4);
            localdof.push_back(local);
        } else {
            owner.push_back(global + addSet[edge]);
            intEdges.push_back(edge);
            localdof.push_back((edge+3)%4 + dofi*4);
        }

    }  // else the dof is supplemental bubble function 
       // The integral domain is then confined to the edge

    // Integrate over domain
    for (int d = 0; d<owner.size(); d++){

         basis_.GetCorners(mi, owner.at(d));
         vertexSet corners = basis_.corners();

         vertexSet corner = {corners.at((intEdges.at(d)+3)%4),
                             corners.at(intEdges.at(d))};

         double len = length(corner);

         for (int g = 0; g<gpe.size(); g++){
             vertex mapped = GaussMapPointsEdge({gpe[g]}, corner);
             vertex evap = br_.ComputeBRmixed(basis_,mapped,localdof.at(d));

             // Evaluate traction on the boundary
             vertex tract = traction(mapped, pp);

             work += len/2.0 * gwe[g] * tract[0]*evap[0] + 
                                        tract[1]*evap[1];  
         }
    }

    return work;
}

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
                         &tmp.val,INSERT_VALUES);
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

PetscErrorCode CreateNeumBndryVec(const int& totalDof,
                                  const int& diriDof,
                                  ReducedSys * resys,
                                  bndryVal& bndryNeum,
                                  bndryVal& bndryDiri){

    // Create Neumann boundary vector in the serial manner

    PetscCall(VecCreate(PETSC_COMM_WORLD, &resys->neum));
    PetscCall(VecSetSizes(resys->neum, PETSC_DECIDE, totalDof - diriDof));
    PetscCall(VecSetUp(resys->neum));

    VecZeroEntries(resys->neum);

    int count = 0;

    // Loop through all degree of freedoms
    for (int globDof=0; globDof<totalDof; globDof++){

        // skip when it is dirichlet dof
        auto keyDiri = bndryDiri.find(globDof);

        if (keyDiri == bndryDiri.end()){

            auto keyFind = bndryNeum.find(globDof);
            if(keyFind != bndryNeum.end()) {
                VecSetValues(resys->neum, 1, &count, &keyFind->second.val, INSERT_VALUES);
            }
            
            count++;
        }
    }

    PetscCall(VecAssemblyBegin(resys->neum));
    PetscCall(VecAssemblyEnd(resys->neum));

    assert(count == totalDof - diriDof);

    return PETSC_SUCCESS;
}
