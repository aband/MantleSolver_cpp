#include "serial_solver.h"

bndryType bMarker(const MeshInfo& mi, std::vector<int> work,
                  const std::string& name,
                  const std::vector<double>& parameter){

    bndryType type = missed;

    if (work[0] != -1){

        if (name == "BDM"){

            // Bndry dof
            indice globalCell = Bend(mi,work[0]); 
            int edge = work[1]%4;
            type = bndryTypeMarkerDarcy(mi, globalCell, edge);

        } else if (name == "BR"){

            indice globalCell = Bend(mi,work[0]);
            type = bndryTypeMarker(mi, globalCell, work[1], parameter);
        }

    }

    return type;
}

double DarcyStokes::AssignBndrySupVal(const vertexSet& edgeCorner,
                                      const vertex& nu,
                                      PhysProperty * pp){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

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

std::array<double, 2> DarcyStokes::AssignBndryValsDarcy(const vertexSet& edgeCorner,
                                                        const vertex& nu,
                                                        const double& len,
                                                        const int& edge,
                                                        PhysProperty * pp){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

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

    work[0] = (d*l0-b*l1)/(a*d-b*b);
    work[1] = (a*l1-b*l0)/(a*d-b*b);

    return work;
}

int DarcyStokes::computeEssenVals(const MeshInfo& mi,
                                  int i, int j, int edge, PhysProperty * pp){

    // Get global cell index
    indice gcell {i,j};

    // Extract corners of this element
    basis_.GetCorners(mi, gcell);
    vertexSet fullCorners = basis_.corners();

    // Extract edge corners and its corresponding normal vector
    vertexSet edgeCorners = {fullCorners.at((edge+3)%4), 
                             fullCorners.at(edge)};

    vertex nu = basis_.unitnormal(edge);

    double len = length(edgeCorners);

    // Stokes ================================================================================

    std::array<int, 12>  elementDOF = br_.LocalToGlobal(mi, gcell);

    vertex essenVal = bndryVs(edgeCorners.at(1), pp);

    double supVal = AssignBndrySupVal(edgeCorners, nu, pp);
//cout << i << "  " << j << "  "  << edge  << " " << supVal << endl << endl;
    std::array<double,3> tmpVal {essenVal[0], essenVal[1], supVal};

    // Store three values
    for (int d=0; d<3; d++){
        int locd = edge + d*4;

        bndryStokesAll.insert(std::make_pair<int, bndryInfo>
                              ((int)elementDOF[locd], {locd, gcell, tmpVal[d], 0}));
    }

    // Darcy  ================================================================================

    std::array<int, 8> elementDOFDarcy = hdiv_.LocalToGlobal(mi, gcell);

    // Get dirichlet boundary value assigned to boundary dofs
    std::array<double, 2> dVals = AssignBndryValsDarcy(edgeCorners, nu, len, edge, pp);
 
    bndryDarcyAll.insert(std::make_pair<int, bndryInfo>
                         ((int)elementDOFDarcy[edge], {edge, gcell, dVals[0], 0})); 

    bndryDarcyAll.insert(std::make_pair<int, bndryInfo>
                         ((int)elementDOFDarcy[edge+4], {edge+4, gcell, dVals[1], 0})); 

    return 1;
}

int DarcyStokes::ComputeEssenBndryAll(const MeshInfo& mi,
                                      PhysProperty * pp,
                                      const std::vector<double>& param){

    // Compute essential boundary condition on every dofs
    // store "right" and supp dof only on each edge

    // left and right edges
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){

       computeEssenVals(mi, 0,j,0, pp); 
       computeEssenVals(mi, mi.MPIglobalCellSize[0]-1,j,2,pp);

    }

    // bottom and top edges
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){
      
       computeEssenVals(mi, i,0,1, pp);
       computeEssenVals(mi, i,mi.MPIglobalCellSize[1]-1,3,pp);

    }

    return 1;
}

int DarcyStokes::computeNaturVals(const MeshInfo& mi, 
                                  int i, int j, int edge, PhysProperty * pp){

    // Get glonal cell index
    indice gcell {i,j};

    // Extract corners of this element
    basis_.GetCorners(mi, gcell);
    vertexSet fullCorners = basis_.corners();

    // Extract edge corners and its corresponding normal vector
    vertexSet edgeCorners = {fullCorners.at((edge+3)%4), 
                             fullCorners.at(edge)};

    vertex nu = basis_.unitnormal(edge);

    double len = length(edgeCorners);

    // Stokes ================================================================================
    // For Stokes part q = p_bar - rho_l g z
    // At phi=0 , q_bar = q_s

    // Evaluate all 12 dofs associated with this perticular cell
    std::array<int, 12>  elementDOFStokes = br_.LocalToGlobal(mi, gcell);

    // Darcy ================================================================================
    // For Darcy part q_l = p_l - rho_l g z
    // Hence it is always zero at natural boundary condition

    std::array<int, 8> elementDOFDarcy = hdiv_.LocalToGlobal(mi, gcell);

    // Integration over edge
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double work = 0.0;
    int globaldof;

    // Only compute normal local dofs
    for (int dof : normal_stokes.at(edge)){

        work = 0.0;

        for (int g=0; g<gwe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},edgeCorners);
 
            std::array<vertex, 12> brval = br_.ComputeBRmixed(basis_, mapped);

            work += len/2.0*gwe[g]*(brval.at(dof)[0]*nu[0]+
                                    brval.at(dof)[1]*nu[1])
                                  * naturvalStokes(mapped, pp);
        }

        // Assign this value to corresponding position
        globaldof = elementDOFStokes.at(dof);
//cout << globaldof << "  " << work << endl;
        bndryStokesAll.at(globaldof).naturval += work;     

    }

    return 1;
}

// Compute Pressure correction on each normal dof
int DarcyStokes::ComputeNaturBndryAll(const MeshInfo& mi,
                                      PhysProperty * pp,
                                      const std::vector<double>& param){

    // Compute essential boundary condition on every dofs
    // store "right" and supp dof only on each edge

    // left and right edges
    for (int j=0; j<mi. MPIglobalCellSize[1]; j++){
    
        computeNaturVals(mi, 0, j, 0, pp);
        computeNaturVals(mi, mi.MPIglobalCellSize[0]-1, j, 2, pp);

    }


    // bottom and top edges
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){
 
        computeNaturVals(mi, i, 0, 1, pp);
        computeNaturVals(mi, i, mi.MPIglobalCellSize[1]-1, 3, pp);

    }

    return 1;
}
