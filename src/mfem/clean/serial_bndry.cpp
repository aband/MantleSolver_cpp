#include "serial_solver.h"

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

int DarcyStokes::ComputeEssenBndryAll(const MeshInfo& mi,
                                      PhysProperty * pp,
                                      const std::vector<double>& param){

    // Compute essential boundary condition on every dofs
    // store "right" and supp dof only on each edge

    // left edge
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
       // Global element index
       int gcell {0,j};

       // Extract corners of this element
       basis_.GetCorners(mi, cell);
       vertexSet fullCorners = basis_.corners();

       // Get global numbering of the dofs associating with this element
       std::array<int, 12> elementDOF = br_.LocalToGlobal(mi, gcell);

       vertexSet edgeCorners = {fullCorners.at(3), 
                                fullCorners.at(0)};

       vertex nu = basis_.unitnormal(0);

             

    }

    // bottom edge


    // right edge


    // top edge

    return 1;
}

int DarcyStokes::ComputeNaturBndryAll(const MeshInfo& mi,
                                      PhysProperty * pp,
                                      const std::vector<double>& param){

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
        // Global element index
        int gcell {0,j};

        // Extract corners of this element
        basis_.GetCorners(mi, cell);
        vertexSet fullCorners = basis_.corners();

        // Get global numbering of the dofs associating with this element
        std::array<int, 12> elementDOF = br_.LocalToGlobal(mi, gcell);

        vertexSet edgeCorners = {fullCorners.at(3), 
                                 fullCorners.at(0)};

         

    }	

    return 1;
}

// Create full list of essential and natural boundary 
// without differentiation of actual boundary type
int DarcyStokes::MarkBndryDOFStokes(const MeshInfo& mi, 
                                    PhysProperty * pp,
                                    const std::vector<double>& param){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    int jstart = mi.MPIlocalCellStart[1];
    int istart = mi.MPIlocalCellStart[0];

    for (int j=jstart; j<jstart + mi.MPIlocalCellSize[1]; j++){
    for (int i=istart; i<istart + mi.MPIlocalCellSize[0]; i++){

        // Global element index
        indice global{i,j};

        // Extract corners of this element
        basis_->GetCorners(mi, global);

        vertexSet fullCorners = basis_->corners();

        // Get global numbering of the dofs associating with this element
        std::array<int, 12> elementDOF = br_->LocalToGlobal(mi, global);

        // Mark all the edges of this element that laying on the boundary
        vector<int> edges;

        markBndryEdge(mi, edges, i, j); 

        for (const auto& edge: edges){
            // Get corners corresponding to this boundary edge
            vertexSet edgeCorners = {fullCorners.at((edge+3)%4),
                                     fullCorners.at(edge)};

            // Get unit normal vector to this boundary edge
            vertex nu = basis_->unitnormal(edge);

            // Get dirichlet boundary nodal value
            // and supplemental bubble function value
            vertex bndryVal = essenbndryVs(edgeCorners[1], pp);

            double supVal = AssignBndrySupVal(edgeCorners, nu, gwe, gpe, pp);

            std::array<double,3> tmpVal {bndryVal[0], bndryVal[1], supVal};

            double neumVal = AssignNaturBndryVal(edgeCorners, nu, gwe, gpe, pp);

        }
    }}

    return 1;
}

int DarcyStokes::MarkBndryDOFDarcy(const MeshInfo& mi,
                                   PhysProperty * pp,
                                    const std::vector<double>& parame){


    return 1;
}
