#include "myFunc.h"
#include "param.h"
#include <random>
#include <petsc.h>

// ================================================================================
inline double InitCD(const vertex& point, PhysProperty * pp){


    if (abs(point[1]) < 120*1000/pp->l0 && abs(point[0]) < abs(point[1]) + pp->l){

        return point[1]*point[1] *0.08 + 0.4;

    } else {
        return 0.2;
    }

}

inline double InitHD(const vertex& point, PhysProperty * pp){

    if (abs(point[1]) < 120*1000/pp->l0 && abs(point[0]) < abs(point[1]) + pp->l){

        // Above Eutectic
        return 1.0;
    } else {
        // Below eutectic
        return -0.05;
    }

}

double ComputePorosity(const vertex& point, Phase * phase){

    double CD = InitCD(point, phase->pp);
    double HD = InitHD(point, phase->pp);

    phase->pPtr->evalPhase(HD, CD);

    return phase->pPtr->phi.mlt;
}

int PorosityOut(double xstart, double ystart, double L, double H, int seed,
                Phase * phase){

    FILE * fp = fopen("InitPoro.dat","w");

    double hx = L/(double)seed;
    double hy = H/(double)seed;

    //phaseState * pPtr = new phaseState();

    for (int j=0; j<seed; j++){
    for (int i=0; i<seed; i++){
        vertex point {xstart + hx*i , ystart + hy*j}; 
        fprintf(fp, "%f ", ComputePorosity(point,phase));
    }fprintf(fp, "\n");}

    fclose(fp);

    return 1;

}

// == Transport boundary settings
// Reconstruction stencils on boundary for WENO reconstructions
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


