// upwelling magma simulation

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

    if (abs(point[1]) < 120*1000/pp->l0 && abs(point[0]) < abs(point[1]) + pp->l){
        double value = 0.05*pow((120*1000/pp->l0 - abs(point[1]))/(120*1000/pp->l0),2) * 
                               (1-abs(point[0])/(abs(point[1])+pp->l));

        return value;
    } else {
        return 0.0;
    }

//    return 0.0;

    // Constant porosity
    //return pp->phi0;
}

double AssignPorosity(double phi_f){
    return 1-phi_f;
}

// ===================================================

const vertex darcyPressureGrad(const vertex& point, PhysProperty * pp){

    // Test case 1: 
    // Balanced pressure (divergence free)
    // return {-1.0,1.0};

    // =========================================================

    // Test case 2:
    // Unbalanced pressure
    // return {2.0,2.0};

    // =========================================================

    // Test case 3:
    // Constant porosity
    return {0.0,0.0}; 
}

const vertex stokesPressureGrad(const vertex& point, PhysProperty * pp){

    // Test case 1:
    // Balanced pressure (divergence free)
    // return {-1.0/pow(pp->phi0,0.5),1.0/pow(pp->phi0,0.5)};

    // =========================================================

    // Test case 2:
    // Unbalanced pressure
    return {0.0,0.0};
}

const vertex divdivVel(const vertex& point, PhysProperty * pp){

    return {2*point[1], -2*point[0]};

}

const vertex stress(const vertex& point, PhysProperty * pp){

    // Calculate deviatoric stress
    double coef = 4*pow(pp->phi0,0.5)/(3*(1-pp->phi0));

    // Test case 1:
    // Balanced pressure
    // return { point[1],
    //         -point[0]};

    // =========================================================

    // Test case 2:
    // Unbalanced pressure
    return {coef, coef};
}

// =============================================================

// Boundary values

vertex bndryVs(const vertex& point, PhysProperty * pp){

    // Stokes
    vertex work {0.0,0.0};
    double coef = 0.0;

    // Test case 1:
    // ====== Test Balanced pressure ======
    //work[0] = point[0]*point[0]*point[1];
    //work[1] = -point[1]*point[1]*point[0];

    // ==================================================

    // Test case 2:
    // ====== Test Unbalanced pressusre ======
    //coef = pow(pp->phi0,0.5)/(1-pp->phi0);
    //work[0] = coef * point[0]*point[0];
    //work[1] = coef * point[1]*point[1];

    // ==================================================

    // Test case 3:
    // Corner Flow
    double x, z;

    if (point[0] < 0.0) {
        x = point[0] - pp->l;
    }else{
        x = point[0] + pp->l;
    }

    z = point[1];

    coef = 2*pp->U0/(3.14159265358979323846*(x*x+z*z))/pp->u0;
    //coef = 2/(3.14159265358979323846*(x*x+z*z));

    work =  {atan(x/z)*(x*x+z*z) - x*z,
             -z*z};

    work *= coef;

    // ==================================================
/*
    // Test Case 4:
    double scale = -1*pp->U0/pp->u0;
    //double scale = 0.002;
    if (point[1] < -0.999*pp->L0/pp->l0){
        work[0] = 0.0;
        work[1] = scale;
    } else if (point[0] < 0 || point[0] > 0){
        work[0] = scale * point[0] / abs(point[0]);
        work[1] = 0.0;
    } else {
        work[0] = 0.0;
        work[1] = 0.0;
    }

    // bottom corner
    //if (point[1] < -0.999/pp->l0 && (point[0] < -0.999/pp->l0 || point[0] > 0.999/pp->l0)){
    //    work[0] = scale * point[0]/abs(point[0]);
    //    work[1] = scale;
    //}
*/

    return work;
}

vertex bndryu(const vertex& point, PhysProperty * pp){

    // Darcy
    vertex work {0.0,0.0};
    double coef = 0.0;

    // Test case 1:
    // ====== Test Balanced pressure ========
    //work[0] = point[0]*point[0]*point[1];
    //work[1] = -point[1]*point[1]*point[0];

    // Test case 2:
    // ====== Test unbalanced pressure ======
    //coef = -1.0/pow(pp->phi0,0.5)/(1-pp->phi0);
    //work[0] = coef * point[0]*point[0];
    //work[1] = coef * point[1]*point[1];


    // Test case 3:
    // Constant porosity
    double x,z;

    if (point[0] < 0.0){
        x = point[0] - pp->l;
    }else {
        x = point[0] + pp->l;
    }

    z = point[1];

    // Point wise porosity
    double phi_f = AssignPorosity(point, pp);

    double rho_r = pp->rho_f*phi_f + pp->rho_s*(1-phi_f);

    //double coef1 = (1-pp->phi0) * pow(pp->phi0,2+2*pp->theta);
    //double coef2 = 4*pp->mu_s*pp->U0/(3.14159265358979323846*(x*x+z*z)*pp->x0*pp->x0) /rho_r /pp->gy;
    double coef1 = (1-phi_f)*pow(phi_f,2+2*pp->theta); 

    double coef2 = 4*pp->U0/pp->u0/(3.14159265358979323846*(x*x+z*z)*(x*x+z*z));

    work[0] = coef1*coef2*2*x*z;
    work[1] = coef1*coef2*(z*z-x*x);

    work[0] += coef1 * 0;
    work[1] += coef1 * 1;

    //cout << point[0] << " " << point[1] << "  " << phi_f << " " << work[1] << endl;

    // Scale the velocity
    //work[0] /= pp->phi0;
    //work[1] /= pp->phi0;

    return work;
}

// ==============================================================
const vertex darcyForce(const vertex& point, PhysProperty * pp){

    // Test case 1,2:
    // Return the arbitrarily defined right hand side
    // source term are determined exactly.
    // return bndryu(point,pp) + pow(pp->phi0,0.5)*darcyPressureGrad(point,pp);

    // ==========================================================
    
    // Test case 3:
    // Constant porosity.
    return {0.0,0.0};
}

const vertex stokesForce(const vertex& point, PhysProperty * pp){

    // Test case 1,2:
    // Source term is determined exactly.
    // return -1*divdivVel(point)+stokesPressureGrad(point);
    // return -1*2*(1-pp->phi0)*stress(point,pp) + stokesPressureGrad(point,pp);
 
    // ==========================================================

    // Test case 3:
    // Constant porosity.
    // Returns nondimensionalized gravity.
    //return {0.0,-1*(1-AssignPorosity(point, pp))/pp->l0};
    return {0.0, -1*(1-AssignPorosity(point, pp))};
    //return {0.0,0.0};
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

    // Current problem setup:
    // inflow dirichlet on the bottom
    // outflow dirichlet on the left and right sides
    // tangential dirichlet on the top side

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

    if (global[1] == mi.MPIglobalCellSize[1]-1 && exit(1.05, mi.MPIglobalCellSize[0], global[0]) ){
        type = neumann; 
    } else {
        type = dirichlet;
    }

    return dirichlet;
}
