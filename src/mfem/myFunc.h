#ifndef MYFUNC_H_
#define MYFUNC_H_

#include "util.h"

typedef struct {

    double theta;
    double mu_s;
    double mu_f;
    double rho_f = 2800;
    double rho_s = 3300;
    double gx = 0.0;
    double gy = -10.0;

    double l = 20;
} PhysProperty;

double AssignPorosity(const vertex& point, const double& l);

// ===================================================
// Define boundary condition
// ===================================================

// True solution
std::array<double, 3> trueSol(const vertex& point);

// Dirichlet value defined on the boundary
const vertex Dirichlet_val(const vertex& point);

// Return source term for darcy system as sum of velocity and pressure gradient  
const vertex darcyForce(const vertex& point);

// Return source term for stokes system as sum of true solutions 
const vertex stokesForce(const vertex& point);

#endif
