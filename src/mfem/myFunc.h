#ifndef MYFUNC_H_
#define MYFUNC_H_

#include "util.h"

typedef struct {

    double theta ;
    double mu_s  ;
    double mu_f  ;
    double rho_f ;
    double rho_s ;
    double gx    ;
    double gy    ;
    double invk0 ;
    double phi0  ;
    double U0    ;

    double x0    ;
    double u0    ;
    double p0    ;

    double l;
} PhysProperty;

double AssignPorosity(const vertex& point, PhysProperty * pp);

void AssignPhyProperties(PhysProperty * pp);
// ===================================================
// Define boundary condition
// ===================================================

// True solution
std::array<double, 3> trueSol(const vertex& point);

// Dirichlet value defined on the boundary
const vertex Dirichlet_val(const vertex& point);

// Return source term for darcy system as sum of velocity and pressure gradient  
const vertex darcyForce(const vertex& point, PhysProperty * pp);

// Return source term for stokes system as sum of true solutions 
const vertex stokesForce(const vertex& point, PhysProperty * pp);

// Boundary Condition
vertex bndryVs(const vertex& point, PhysProperty * pp);

vertex bndryu(const vertex& point, PhysProperty * pp);

#endif
