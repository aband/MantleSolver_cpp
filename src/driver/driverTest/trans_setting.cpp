#include "trans_param.h"
#include "coupled.h"

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
       (globalCell[0] == 0 && globalCell[1] == 0) ||
       (globalCell[0] == mi.MPIglobalCellSize[0]-1 && globalCell[1] == mi.MPIglobalCellSize[1]-1) 
      ){
       return true;
   } else {
       return false;
   }
}

// Boundary categary

std::string WEAK_COUPLED::transport::location(const MeshInfo& mi,
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

double InitCD(const valarray<double>& point,
              const vector<double>& param){

    return 0;
}

double InitHD(const valarray<double>& point,
              const vector<double>& param){

    return 0;
}
