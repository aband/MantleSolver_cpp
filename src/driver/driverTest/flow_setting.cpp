#include "myFunc.h"

void AssignPhyProperties(PhysProperty * pp){

    pp->theta = 0.0;
    pp->mu_s  = 1e19;
    pp->mu_f  = 1.0;
    pp->rho_f = 3000;
    pp->rho_s = 3000;
    pp->gx    = 0.0;
    pp->gy    = -10.0;
    pp->invk0 = 1.0/(1e-8);
    pp->phi0  = 0.4;
    pp->U0    = 1e-9;
    pp->L0    = 160*1000;

    // Non dimensionalization parameters

    double rho_r = pp->rho_f*pp->phi0 + 
                   pp->rho_s*(1-pp->phi0); 

    pp->l0    = pow(pp->mu_s/pp->invk0/pp->mu_f,0.5);
    pp->p0    = pp->gy*pp->l0*rho_r;
    pp->u0    = pp->gy*rho_r/pp->mu_f/pp->invk0;

    pp->l = 20/pp->l0;
}

double AssignPorosity(const vertex& point, PhysProperty * pp){

    // used to identify incorrect porosity
    return -1.0;
}

double AssignPorosity(double phi_f){
    return 1-phi_f;
}

// ===================================================
// Below not needed if no true solution posted
/*
const vertex darcyPressureGrad(const vertex& point, PhysProperty * pp){

    return {0.0,0.0}; 
}

const vertex stokesPressureGrad(const vertex& point, PhysProperty * pp){

    return {0.0,0.0};
}

const vertex divdivVel(const vertex& point, PhysProperty * pp){

    return {2*point[1], -2*point[0]};
}

const vertex stress(const vertex& point, PhysProperty * pp){

    // Calculate deviatoric stress
    double coef = 4*pow(pp->phi0,0.5)/(3*(1-pp->phi0));
}
*/
// =============================================================

// Boundary values
// Constant upwelling velocity ascending model
vertex bndryVs(const vertex& point, PhysProperty * pp){

    // Stokes
    double V0 = 3.2/100/(365*24*60*60); //3.2 (cm/y)

    V0 = V0 / pp->u0 * -1;

    return {0.0,0.01};
}

vertex bndryu(const vertex& point, PhysProperty * pp){

    // Darcy

    return {0.0,0.0};
}

// ==============================================================
const vertex darcyForce(const vertex& point, PhysProperty * pp){

    // Constant porosity.
    return {0.0,0.0};
}

const vertex stokesForce(const vertex& point, PhysProperty * pp){

    // Returns nondimensionalized gravity.
    // Attention!!! It should not be scaled by porosity
	 // porosity scale will be added in another function
    return {0.0, -0.0};
}

const vertex traction(const vertex& point, PhysProperty * pp){
    // return traction defined on the boundary
/*
    if (point[0] < -0.999*pp->L0/pp->l0 || point[0] > 0.999*pp->L0/pp->l0) {
        return {0.0,abs(point[1])*pp->L0/pp->l0};
    } else {
        return {0.0,0.0};
    }
*/
    return {0.0,0.0}; // free stress
}

// =========================================================================
const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global,
                                const int& local){

    bndryType type = missed;

    set<int> top_normal {11,7,6};
    set<int> left_normal {0,3,8};
    set<int> right_normal {1,2,10};
    set<int> bottom_normal {4,5,9};

/*
    if (global[1] == mi.MPIglobalCellSize[1] - 1 && global[0] != 0){
        if (local == 11 || local == 7){
            type = neumann;
        } else {
            type = dirichlet;
        }

    }else if (global[1] == mi.MPIglobalCellSize[1] - 1 && global[0] == 0){

        if (local == 11){
            type = neumann;
        } else {
            type = dirichlet;
        }

    } else {
        type = dirichlet;
    }
*/

    if (global[1] == 0 ){

        type = dirichlet;

    } else if (global[0] == 0){

        if (local == ){

        }

    } else if (global[0] == mi.MPIglobalCellSize[0] - 1){

    } else {

        type = neumann;
    }

    return type;
}

bool exit(double range, int M, int i){

    if (i < M/2 + range && i > M/2 - range){
        return true;
    } else {
        return false;
    }

}

const bndryType bndryTypeMarkerDarcy(const MeshInfo& mi,
                                     const indice& global,
                                     const int& edge){

    bndryType type = missed;

    //if (global[1] == mi.MPIglobalCellSize[1]-1 && exit(1.05, mi.MPIglobalCellSize[0], global[0]) ){
    if (global[1] == mi.MPIglobalCellSize[1]-1){
        if (edge == 3){
            type = neumann; 
        }else {
            type = dirichlet;
        }
    } else {
        type = dirichlet;
    }

    return type;
}
