#include "myFunc.h"

// The input file for corner flow problem
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

    double rho_r = 3000; 

    // Non dimensionalization parameters
    pp->l0    = pow(pp->mu_s/pp->invk0/pp->mu_f,0.5);
    pp->p0    = pp->gy*pp->l0*rho_r;
    pp->u0    = pp->gy*rho_r/pp->mu_f/pp->invk0;

    pp->l = 20/pp->l0;
}


