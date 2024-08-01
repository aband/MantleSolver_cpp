#include "myFunc.h"
#include "param.h"
#include <petsc.h>

// ================================================================================
inline double InitCD(const vertex& point, PhysProperty * pp){

    if (abs(point[1]) < 120*1000/pp->l0 && abs(point[0]) < abs(point[1]) + pp->l){

        return 0.4;
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

    int state = phase->pPtr->EvalPhaseRegion(CD,HD);
    phase->pPtr->EvalPhase(state, CD, HD);

    return phase->pPtr->phi.fluid;
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
