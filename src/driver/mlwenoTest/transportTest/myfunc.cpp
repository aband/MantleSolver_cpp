#include "myfunc.h"

/**
 * Change functions for transport part here.
 * Define transport functions and derivatives
 * 2D Burger's equation
 *
 */
const std::array<double, 2> advFunc(double u){
    return {0.5*u*u,0.5*u*u}; 
}

const std::array<double, 2> dAdvFunc(double u){
    return {u,u};
}

// ======== Diffusion =======================
double diffFunc(double u){
    return u;
}

double dDiffFunc(double u){
    return 1;
}

// ======== Initial distribution ============
double Distribution(const vertex& point, 
                    const vector<double>& param){

//    if (point[0] < 0){
//        return -1;
//    } else {
//        return 1;
//    }

	 //if (point[0]<-1.0/param[0]){
//		  return point[0]*point[0]+point[1]*point[1];
//	     return sin(point[0]*3.0)+cos(point[1]/2.0) + point[0]*(point[1]+1);
//	 } else {
//		  return point[0]*point[0]*point[1]*point[1] + 1.0;
//	     return sin(point[0]*3.0)+cos(point[1]/2.0) + point[0]*(point[1]+1) + 1;
//	 }

    //return sin(point[0]*3.0+0.5)+cos(point[1]/2.0-0.2) + pow(point[0]+0.1,3)*(point[1]+1);
    //return point[0]*point[0] + point[1]*point[1];
    //return point[0] + point[1];

    // Initial value for sine wave 2D Burger's equation
    return pow(sin(M_PI*(point[0]+1)/2),2)*pow(sin(M_PI*(point[1]+1)/2),2);

    //if (abs(point[0])+abs(point[1])<0.5){
    //    return 1;
    //} else {
    //    return 0;
    //}
}

// ======== Split different region ============

bool left_boundary(const indice& globalCell,
                   const MeshInfo& mi){

    if (globalCell[0] == 0 &&
        globalCell[1] > 0 && 
        globalCell[1] < mi.MPIglobalCellSize[1]-1) {
        return true;
    } else {
        return false;
    }
}

bool right_boundary(const indice& globalCell,
                    const MeshInfo& mi){

    if (globalCell[0] == mi.MPIglobalCellSize[0]-1 &&
        globalCell[1] > 0 &&
        globalCell[1] < mi.MPIglobalCellSize[1]-1){
        return true;
    } else {
        return false;
    }

}

bool top_boundary(const indice& globalCell,
                  const MeshInfo& mi){

    if (globalCell[1] == mi.MPIglobalCellSize[1]-1 &&
        globalCell[0] > 0&&
        globalCell[0] < mi.MPIglobalCellSize[0]-1){
        return true;
    } else {
        return false;
    }

}

bool bottom_boundary(const indice& globalCell,
                     const MeshInfo& mi){

    if (globalCell[1] == mi.MPIglobalCellSize[1]-1 &&
        globalCell[0] > 0 &&
        globalCell[0] < mi.MPIglobalCellSize[0]-1){
        return true;
    } else {
        return false;
    }

}

bool interior(const indice& globalCell, 
              const MeshInfo& mi){
   if (globalCell[0] > 0 && globalCell[0] < mi.MPIglobalCellSize[0]-1 &&
       globalCell[1] > 0 && globalCell[1] < mi.MPIglobalCellSize[1]-1){
       return true;
   } else {
       return false;
   }
}

Location assignLocation(const indice& globalCell,
                        const MeshInfo& mi){

   if (left_boundary(globalCell, mi)){
       return leftBndry;
   } else if (right_boundary(globalCell, mi)){
       return rightBndry;
   } else if (top_boundary(globalCell, mi)){
       return topBndry;
   } else if (bottom_boundary(globalCell, mi)){
       return bottomBndry
   } else {
       return interior;
   }

}
