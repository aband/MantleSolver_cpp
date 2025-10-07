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

bool top(const indice& globalCell,
         const MeshInfo& mi){


    if (globalCell[1] > mi.MPIglobalCellSize[1]-1 - 4){
        return true;
    } else {
        return false;
    }

}

// Boundary categary

std::string location(const MeshInfo& mi,
                     const indice& globalCell){

    std::string loc = "wrong";
/*
    if (interior(globalCell, mi)){

       loc = "interior";

    } else if (edge(globalCell, mi)){

       loc = "edge";

    } else if(corner(globalCell, mi)){

       loc = "corner";

    }
*/

    if (top(globalCell, mi)){
       loc = "top";
    } else {
       loc = "all";
    }
    return loc;
    //return "all";
}

// Initialize dimensionless composition and enthalpy
double InitCD(const valarray<double>& point,
              const vector<double>& param){

    // Constant composition value 
    return 0.15;
}

double InitHD(const valarray<double>& point,
              const vector<double>& param){

    // Linear simple distribution of enthalpy
	 // We pass nondimensionalize normalization factor in param.at(0)
    //double HD = 0.01;

    double cut  = -0.2;
    cut = -1.5;
    double slop = 2.8;
    double head = 2.9;

    double HD = head - slop*cut -slop*pow(point[1]-cut,1);

    if (point[1] < cut){HD = head - slop*cut;}

    return HD;
}

double advfunc(const double& u, 
               const vertex& vel, const vertex& unitnormal){

//    cout << "Called correct one" << endl;

    return u *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
}

double dfdu(const double& u){

    return 1.0;
}

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work){

    // compute df/du = df/dR * dR/du

    double direction = dfdu(u)*(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);

    return 1;
}
