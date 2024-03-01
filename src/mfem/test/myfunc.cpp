#include "myFunc.h"

double AssignPorosity(const vertex& point, const double& l){

    //if (point[1] < 120000 && abs(point[0]) < point[1] + l){
    //    return 0.05*pow(1.0-point[1]/120000,2) * (1-abs(point[0])/(l+point[1]));
    //} else {
    //    return 0.0;
    //}

    return 0.5;
}

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

    // Non dimensionalization parameters

    double rho_r = pp->rho_f*pp->phi0 + 
                   pp->rho_s*(1-pp->phi0); 

    pp->x0    = pow(pp->mu_s/pp->invk0/pp->mu_f,0.5);
    pp->p0    = pp->gy*pp->x0*rho_r;
    pp->u0    = pp->gy*rho_r/pp->mu_f/pp->invk0;

    pp->l = 2.0;
}

// ===================================================

const vertex darcyPressureGrad(const vertex& point, PhysProperty * pp){

    // Test balanced pressure
    //return {-1.0,1.0};

    // Test unbalanced pressure
    return {0.0,0.0};
}

const vertex stokesPressureGrad(const vertex& point, PhysProperty * pp){

    // Test balanced pressure
    //return {-1.0/pow(pp->phi0,0.5),1.0/pow(pp->phi0,0.5)};

    // Test unbalanced pressure
    return {2,2};

}

const vertex divdivVel(const vertex& point, PhysProperty * pp){

    return {2*point[1], -2*point[0]};

}

const vertex stress(const vertex& point, PhysProperty * pp){

    // Calculate deviatoric stress
    double coef = 4*pow(pp->phi0,0.5)/(3*(1-pp->phi0));

    // Test balanced pressure
    //return { point[1],
    //        -point[0]};

    // Test unbalanced pressure
    return {coef, coef};
}

// ====================================================================

// Boundary values

vertex bndryVs(const vertex& point, PhysProperty * pp){
    // Stokes

    // Rewrite it with non dimensionalized versioin

    vertex work {0.0,0.0};

    double x, z;

    if (point[0] < 0.0) {
        x = point[0] - pp->l/pp->x0;
        x*= -1;
    }else{
        x = point[0] + pp->l/pp->x0;
    }

    z = point[1];

    double coef = 2*pp->U0/(3.14159265358979323846*(x*x+z*z))/pp->u0;

    work =  {atan(x/z)*(x*x+z*z) - x*z,
             -z*z};

    work *= coef;

    // ====== Test Balanced pressure ======
    //work[0] = point[0]*point[0]*point[1];
    //work[1] = -point[1]*point[1]*point[0];
 
    // ====== Test Unbalanced pressusre ======
    coef = pow(pp->phi0,0.5)/(1-pp->phi0);
    work[0] = coef * point[0]*point[0];
    work[1] = coef * point[1]*point[1];

    return work;
}

vertex bndryu(const vertex& point, PhysProperty * pp){
    // Darcy

    vertex work {0.0,0.0};

    double x,z;

    if (point[0] < 0.0){
        x = point[0] - pp->l/pp->x0;
        x *= -1;
    }else {
        x = point[0] + pp->l/pp->x0;
    }

    z = point[1];

    double rho_r = pp->rho_f*pp->phi0 + pp->rho_s*pp->phi0;

    double coef1 = (1-pp->phi0) * pow(pp->phi0,2+2*pp->theta);
    double coef2 = 4*pp->mu_s*pp->U0/(3.14159265358979323846*(x*x+z*z)*pp->x0*pp->x0) /rho_r /pp->gy;

    work[0] = coef1*coef2*2*x*z;
    work[1] = coef1*coef2*(pow(z,2) - pow(x,2));

    work[0] += coef1 * 0;
    work[1] += coef1 * 1;

    // ====== Test Balanced pressure ========
    //work[0] = point[0]*point[0]*point[1];
    //work[1] = -point[1]*point[1]*point[0];

    // ====== Test unbalanced pressure ======
    double coef = -1.0/pow(pp->phi0,0.5)/(1-pp->phi0);
    work[0] = coef * point[0]*point[0];
    work[1] = coef * point[1]*point[1];

    return work;
}

// ==============================================================
const vertex darcyForce(const vertex& point, PhysProperty * pp){

    // Return the arbitrarily defined right hand side
    // source term.
    return bndryu(point,pp) + pow(pp->phi0,0.5)*darcyPressureGrad(point,pp);
}

const vertex stokesForce(const vertex& point, PhysProperty * pp){

    //return -1*divdivVel(point)+stokesPressureGrad(point);
    return -1*2*(1-pp->phi0)*stress(point,pp) + stokesPressureGrad(point,pp);
}


