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
    pp->phi0  = 0.5;
    pp->U0    = 1e-9;
    pp->x0    = 160*1000;

    pp->l = 20.0;
}

// ===================================================

std::array<double, 3> trueSol(const vertex& point){

    array<double, 3> work;

    // ======================================================
    // Darcy test problem
    // return a predefined true solution
    // return <ux, uy, p> in this order
    // A manufactured solution satisfying Darcy equation

    // First scenerio
    // u + grad(p) = 0
    // div(u)      = 0
    // ux = -x/(x^2+y^2)
    // uy = -y/(x^2+y^2)
    // p  = 1/2 ln(x^2+y^2)

    //work[0] = -point[0]/(point[0]*point[0] + point[1]*point[1]);
    //work[1] = -point[1]/(point[0]*point[0] + point[1]*point[1]);
    //work[2] = 0.5*log(point[0]*point[0] + point[1]*point[1]);

    // Second scenerio
    // Divergence free linear velocity with arbitrary defined pressure field

    // Third scenerio
    work[0] = pow(point[0],2)*point[1];
    work[1] = -pow(point[1],2)*point[0];
    work[2] = -point[0] + point[1];

    // =================================================================
    // Test for Stokes problem

    // Constant true solution
    //work[0] = pow(point[0],3)*pow(point[1],2);
    //work[1] = -pow(point[1],3)*pow(point[0],2);
    //work[0] = point[1]*point[1];
    //work[0] = pow(point[1],2);
    //work[1] = 0;
    //work[2] = 0.0;

    return work;
}

const vertex darcyPressureGrad(const vertex& point){

    // Auxiliary function.
    // Returns the gradient of scalar pressure field
    return {0.0,0.0};
}

const vertex stokesPressureGrad(const vertex& point){

    //return {cos(point[0])*sin(point[1]),
    //        sin(point[0])*cos(point[1])};

    return {-1.0,1.0};
    //return {-point[0], point[1]};
}

const vertex divdivVel(const vertex& point){

    //return {-2*cos(point[0])*sin(point[1]),
    //         2*sin(point[0])*sin(point[1])};

    return {6*point[0]*pow(point[1],2) + 2*pow(point[0],3),
            -6*point[1]*pow(point[0],2) - 2*pow(point[1],3)};
    //return {2,0.0};
    //return {2*point[1], -2*point[0]};

}

const vertex Dirichlet_val(const vertex& point){

    std::array<double, 3> truesol = trueSol(point);

    // Test of Darcy part
    return {truesol[0], truesol[1]};
}

const vertex darcyForce(const vertex& point){

    // Return the arbitrarily defined right hand side
    // source term.
    std::array<double, 3> truesol = trueSol(point);

    vertex gradpressure = darcyPressureGrad(point);

    //return {truesol[0] + gradpressure[0], truesol[1] + gradpressure[1]};

    return {0.0,0.0};
}

const vertex stokesForce(const vertex& point){

    return -1*divdivVel(point)+stokesPressureGrad(point);
}

// Boundary values

vertex bndryVs(const vertex& point, PhysProperty * pp){

    // Rewrite it with non dimensionalized versioin

    vertex work {0.0,0.0};

    double x, z;

    if (point[0] < 0.0) {
        x = point[0] - pp->l;
    }else{
        x = point[0] + pp->l;
    }

    z = point[1];

    double coef = 2/(3.14159265358979323846*(x*x+z*z));

    work =  {atan(x/z)*(x*x+z*z) - x*z,
             -z*z};

    work *= coef;

    return work;
}

vertex bndryu(const vertex& point, PhysProperty * pp){

    vertex work {0.0,0.0};

    double coef1 = 1.0/pp->invk0 * (1-pp->phi0) * pow(pp->phi0,2+2*pp->theta)/ pp->mu_f;
    double coef2 = 4*pp->mu_s/(3.14159265358979323846*(point[0]*point[0]+point[1]*point[1]));

    double rho_r = pp->rho_f*pp->phi0 + pp->rho_s*pp->phi0;

    work[0] = coef1*coef2*2*point[0]*point[1]/pp->x0/pp->x0;
    work[1] = coef1*coef2*(pow(point[1],2) - pow(point[0],2))/pp->x0/pp->x0;

    work[0] += coef1 * rho_r * pp->gx;
    work[0] += coef1 * rho_r * pp->gy;

    return work;
}
