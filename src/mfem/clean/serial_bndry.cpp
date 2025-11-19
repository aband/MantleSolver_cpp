#include "serial_solver.h"

static int markBndryEdge(const MeshInfo& mi, 
                         vector<int>& edges,
                         int i, int j){

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

    return 1;
}

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
        basis_.GetCorners(mi, global);

        vertexSet fullCorners = basis_.corners();

        // Get global numbering of the dofs associating with this element
        std::array<int, 12> elementDOF = br_.LocalToGlobal(mi, global);

        // Mark all the edges of this element that laying on the boundary
        vector<int> edges;

        markBndryEdge(mi, edges, i, j); 

        for (const auto& edge: edges){
            // Get corners corresponding to this boundary edge


        }

    }}

    return 1;
}

int DarcyStokes::MarkBndryDOFDarcy(const MeshInfo& mi,
                                   PhysProperty * pp,
                                    const std::vector<double>& parame){


    return 1;
}
