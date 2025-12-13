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
    pp->V0    = 3.2/100/(365*24*60*60); //3.2 (cm/y)

    // Non dimensionalization parameters

    pp->rho_r = pp->rho_s - pp->rho_f;

    pp->l0    = pow(pp->mu_s/pp->invk0/pp->mu_f,0.5);
    pp->p0    = pp->gy*pp->l0*pp->rho_r;
    pp->u0    = pp->gy*pp->rho_r/pp->mu_f/pp->invk0;

    pp->l = 20/pp->l0;
}

// Essential boundary Boundary values =================================
// Constant upwelling velocity ascending model
vertex bndryVs(const vertex& point, PhysProperty * pp){

    // Stokes
    double V0 = pp->U0 / pp->u0 * -1;

    return {0.0,V0};
    //return {0.0, 1.0};
}

vertex bndryu(const vertex& point, PhysProperty * pp){

    // Darcy

    //return {abs(point[1]), 0.0};
    return {0.0,0.0};
}

// ====================================================================
const vertex darcyForce(const vertex& point, PhysProperty * pp){

    // Constant porosity.
    return {0.0,0.0};
}

const vertex stokesForce(const vertex& point, PhysProperty * pp){

    // Returns nondimensionalized gravity.
    // Attention!!! It should not be scaled by porosity
	 // porosity scale will be added in another function
double V0 = pp->V0 / pp->u0;	
    //return {0.0, -1.0/V0};
    return {0.0,-1.0};
}

// ====================================================================
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


