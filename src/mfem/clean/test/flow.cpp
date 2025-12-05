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

    //double rho_r = pp->rho_f*pp->phi0 + 
    //               pp->rho_s*(1-pp->phi0); 

    //double rho_r = pp->rho_s - pp->rho_f;

    double rho_r = pp->rho_s;

    pp->l0    = pow(pp->mu_s/pp->invk0/pp->mu_f,0.5);
    pp->p0    = pp->gy*pp->l0*rho_r;
    pp->u0    = pp->gy*rho_r/pp->mu_f/pp->invk0;

    pp->l = 20/pp->l0;
}

// Essential boundary Boundary values
// Constant upwelling velocity ascending model
vertex bndryVs(const vertex& point, PhysProperty * pp){

    // Stokes
    double V0 = pp->U0 / pp->u0 * -1;

    return {0.0,V0};
    //return {0.0, 1.0};
}

vertex bndryu(const vertex& point, PhysProperty * pp){

    // Darcy

    return {abs(point[1]), 0.0};
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
