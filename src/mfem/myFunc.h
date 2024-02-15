#ifndef MYFUNC_H_
#define MYFUNC_H_

#include "util.h"

typedef struct {

    double theta = 0.0;
    double mu_s  = 10e19;
    double mu_f  = 1.0;
    double rho_f = 2800;
    double rho_s = 3300;
    double gx    = 0.0;
    double gy    = -10.0;
    double invk0 = 1.0/(10e-8);
    double phi0  = 0.5;
    double U0    = 10e-9;

    double l = 20.0;
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

// Boundary Condition
std::array<double, 2> bndryVs(const vertex& point, PhysProperty * pp);

std::array<double, 2> bndryu(const vertex& point, PhysProperty * pp);

#endif
