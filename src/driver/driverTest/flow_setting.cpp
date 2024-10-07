#include "myFunc.h"

void AssignPhyProperties(PhysProperty * pp){

    pp->theta = 0.0;
    pp->mu_s  = 1e19;
    pp->mu_f  = 1.0;
    pp->rho_f = 2800;
    pp->rho_s = 3300;
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
    return {0.0,1.0};
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
    return {0.0, -1*(1-AssignPorosity(point, pp))};
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
/*
const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global, 
                                const int& local){
   
    // Mark boundary condition on each Stokes DOFs
    bndryType type = missed;

    if (global[0] == 0){
        // left side
        if (local == 0 || local == 8){
            type = dirichlet;
        } else if (local == 4){
            type = neumann;
        } 
    } 

    if (global[0] == mi.MPIglobalCellSize[0]-1){
        // right side
        if (local == 2 || local == 10){
            type = dirichlet;
        } else if (local == 6){
            type = neumann;
        } 
    } 

    if (global[1] == 0){
        // bottom side
        if (local == 5 || local == 9){
            type = dirichlet;
        } else if (local == 1){
            type = neumann;
        }
        type = dirichlet;
    }
   
    if (global[1] == mi.MPIglobalCellSize[1]-1){
        // top side
        if (local == 3){
            type = dirichlet;
        } else if (local == 11 || local == 7){
            type = neumann;
        }
    }

    // Bottom two dofs are dealt with separately
    if (global[0] == 0 && global[1] == 0){
        // bottom left
        if (local == 0 || local == 4){
            type = dirichlet;
        }
    }

    if (global[0] == mi.MPIglobalCellSize[1]-1 && global[1] == 0){
        if (local == 1 || local == 5){
            type = dirichlet;
        }
    }

    // the dof is missed during the marking process
    return dirichlet;
} */

const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global,
                                const int& local){

    bndryType type = missed;

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

    return dirichlet;
}

bool exit(double range, int M, int i){

    if (i < M/2 + range && i > M/2 - range){
        return true;
    } else {
        return false;
    }

}

const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global){

    bndryType type = missed;

    //if (global[1] == mi.MPIglobalCellSize[1]-1 && exit(1.05, mi.MPIglobalCellSize[0], global[0]) ){
    if (global[1] == mi.MPIglobalCellSize[1]-1){
        type = neumann; 
    } else {
        type = dirichlet;
    }

    return type;
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

    return dirichlet;
}


