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

// Used for testing two separate darcy and stokes problems
double AssignPorosity(const vertex& point, PhysProperty * pp){

    return 1.0;
}

double AssignPorosity(double phi_f){

    return 1.0;
}

const vertex darcyPressureGrad(const vertex& point, PhysProperty * pp){

    // Test case 1: 
    // Balanced pressure (divergence free)
     return {-1.0,1.0};
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

    return { 3*point[0]*point[1]*point[1] + point[0]*point[0]*point[0],
            -3*point[1]*point[0]*point[0] - point[1]*point[1]*point[1]};
}

const vertex stokesPressureGrad(const vertex& point, PhysProperty * pp){

    // Test case 1:
    // Balanced pressure (divergence free)
    // return {-1.0/pow(pp->phi0,0.5),1.0/pow(pp->phi0,0.5)};
    return {-1.0, 1.0};
}

// Boundary values
vertex bndryVs(const vertex& point, PhysProperty * pp){

	 // Stokes
	 vertex work {0.0,0.0};
	 double coef = 0.0;

	 // Test case 1:
	 // ====== Test Balanced pressure ======
	 work[0] = point[0]*point[0]*point[0]*point[1]*point[1];
	 work[1] = -point[1]*point[1]*point[1]*point[0]*point[0];

    return work;
}

vertex bndryu(const vertex& point, PhysProperty * pp){

    // Darcy
    vertex work {0.0,0.0};
    double coef = 0.0;

    // Test case 1:
    // ====== Test Balanced pressure ========
    work[0] = point[0]*point[0]*point[1];
    work[1] = -point[1]*point[1]*point[0];

    return work;
}

const vertex darcyForce(const vertex& point, PhysProperty * pp){

    // Test case 1,2:
    // Return the arbitrarily defined right hand side
    // source term are determined exactly.
     return bndryu(point,pp) + darcyPressureGrad(point,pp);
}

const vertex stokesForce(const vertex& point, PhysProperty * pp){

    // Test case 1,2:
    // Source term is determined exactly.
    // return -1*divdivVel(point)+stokesPressureGrad(point);
    // return -1*2*(1-pp->phi0)*stress(point,pp) + stokesPressureGrad(point,pp);
    return -2*stress(point,pp) + stokesPressureGrad(point,pp);

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

// For this stokes test problem,
// we test full dirichlet problem.
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
