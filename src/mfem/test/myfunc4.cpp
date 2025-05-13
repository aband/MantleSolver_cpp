#include "myFunc.h"

// Degenerate porosity Darcy problem

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

    //double r = pow(point[0]*point[0] + point[1]*point[1],0.5);
    //double PI = 3.14159265358979323846;

    //if ( r < PI/(2.0*5.0)){
    //    return pow(cos(5.0*r),4);
    //} else {
    //    return 0.0;
    //}

    if (point[0] > 0){
        return 0.4*pow(point[0],4);
    } else {
        return 0.0; 
    }

}


double AssignPorosity(double phi_f){

    return 1.0 - phi_f;
}

// Stokes exact values
vertex bndryVs(const vertex& point, PhysProperty * pp){

    vertex work {0.0,0.0};

    work[0] =  point[0]*point[0]*point[1];
    work[1] = -point[1]*point[1]*point[0];

    return work;
}  

// Darcy exact values
vertex bndryu(const vertex& point, PhysProperty * pp){

    if (point[0] > 0){
        return {point[0]*point[0]/6, -point[0]*point[1]}; 
    } else {
        return {0.0,0.0};
    }
}

// Source term
const vertex darcyForce(const vertex& point, PhysProperty * pp){

    return bndryu(point, pp);
}

const vertex stokesForce(const vertex& point, PhysProperty * pp){

    return {0.0,0.0};

}

const vertex traction(const vertex& point, PhysProperty * PP){

    return {0.0,0.0};
}

const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global, 
                                const int& local){

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

    return dirichlet;
}

const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global,
                                const int& local,
                                const std::vector<double>& parameter){

    return dirichlet;
}

const bndryType bndryTypeMarkerDarcy(const MeshInfo& mi,
                                     const indice& global,
                                     const int& edge){

    return dirichlet;
}
