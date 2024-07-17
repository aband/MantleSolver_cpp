#ifndef SIMPLEPHASE_H_
#define SIMPLEPHASE_H_

// A simplified phase code
#include <math.h>
#include <iostream>

typedef struct{

    double L;  // Latent heat
    double Xe; // 
    double T1; // Melting temperature of phase 1
    double Te; // Eutectic temperature

    double DT; // Temperature difference of T1 and Te

    // Density
    double rho_ice;
    double rho_sal;
    double rho_bri;

    double rho_bi;
    double rho_si;

    // Heat capacity
    double cp_ice;
    double cp_sal;
    double cp_bri;

    double cp_bi;
    double cp_si;

    double Ste;

    // Thermal conductivity
    double kappa_ice;
    double kappa_sal;
    double kappa_bri;

    double kappa_ii;
    double kappa_bi;
    double kappa_si;

    double sysD;

} Phase;

typedef struct {

    double ice;
    double sal;
    double bri;

} hD;

typedef struct {

    double ice;
    double sal;
    double bri;
 
} Phi;

#endif
