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

    return 1;
}

int DarcyStokes::MarkBndryDOFDarcy(const MeshInfo& mi,
                                   PhysProperty * pp,
                                    const std::vector<double>& parame){


    return 1;
}
