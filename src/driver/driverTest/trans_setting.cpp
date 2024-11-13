#include "trans_param.h"

// Location functions
bool interior(const indice& globalCell, 
              const MeshInfo& mi){
   if (globalCell[0] > 0 && globalCell[0] < mi.MPIglobalCellSize[0]-1 &&
       globalCell[1] > 0 && globalCell[1] < mi.MPIglobalCellSize[1]-1){
       return true;
   } else {
       return false;
   }
}

bool edge(const indice& globalCell,
          const MeshInfo& mi){

   // return four edges
   if (// left edge
       (globalCell[0] == 0 && 
        globalCell[1] != 0 && 
        globalCell[1] != mi.MPIglobalCellSize[1]-1) ||
       // right edge
       (globalCell[0] == mi.MPIglobalCellSize[0]-1 && 
        globalCell[1] != 0 && 
        globalCell[1] != mi.MPIglobalCellSize[1]-1) ||
       // bottom edge
       (globalCell[1] == 0 && 
        globalCell[0] != 0 && 
        globalCell[0] != mi.MPIglobalCellSize[0]-1) ||
       // top edge
       (globalCell[1] == mi.MPIglobalCellSize[1]-1 && 
        globalCell[0] != 0 && 
        globalCell[0] != mi.MPIglobalCellSize[0]-1) 
      ){
       return true;
   } else {
       return false;
   }

}

bool corner(const indice& globalCell, 
            const MeshInfo& mi){

   // return four corners

   if ((globalCell[0] == 0 && globalCell[1] == 0) ||
       (globalCell[0] == 0 && globalCell[1] == mi.MPIglobalCellSize[1]-1) ||
       (globalCell[0] == mi.MPIglobalCellSize[0]-1 && globalCell[1] == 0) ||
       (globalCell[0] == mi.MPIglobalCellSize[0]-1 && globalCell[1] == mi.MPIglobalCellSize[1]-1) 
      ){
       return true;
   } else {
       return false;
   }
}

// Boundary categary

std::string location(const MeshInfo& mi,
                     const indice& globalCell){

    std::string loc = "wrong";

    if (interior(globalCell, mi)){

       loc = "interior";

    } else if (edge(globalCell, mi)){

       loc = "edge";

    } else if(corner(globalCell, mi)){

       loc = "corner";

    }

    return loc;
}

// Initialize dimensionless composition and enthalpy
double InitCD(const valarray<double>& point,
              const vector<double>& param){

    // Constant composition value 
    return 0.1;
}

double InitHD(const valarray<double>& point,
              const vector<double>& param){

    // Linear simple distribution of enthalpy
	 // We pass nondimensionalize normalization factor in param.at(0)
    double HD = 0.15;

    HD -= 0.0000005*point[1]*param.at(0);

    return HD;
}

bndryTypeTrans AssignBoundary(const indice& globalCell, const int& edgetype, const std::string& field){

    return noFlow;;
}
