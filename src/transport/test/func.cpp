#include "func.h"

/**
 * Change functions for transport part here.
 * Define transport functions and derivatives
 * 2D Burger's equation
 *
 */
double funcX(vertex x, double u, double t){
    return 0.5*u*u;
}

double dfuncX(vertex x, double u, double t){
    return u;
}

double funcY(vertex x, double u, double t){
    return 0.5*u*u;
}

double dfuncY(vertex x, double u, double t){
    return u;
}

// ======== Diffusion =======================
double diffFunc(double u){
    return u;
}

double dDiffFunc(double u){
    return 1;
}
