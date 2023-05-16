#include "func.h"

/**
 * Change functions for transport part here.
 * Define transport functions and derivatives
 * 2D Burger's equation
 *
 */
double funcX(vertex x, double u, double t){
    return u*u/(u*u+(1-u)*(1-u));
}

double dfuncX(vertex x, double u, double t){
    return (2*u*(1-u))/pow(u*u+(1-u)*(1-u),2);
}

double funcY(vertex x, double u, double t){
    return funcX(x,u,t) *(1-5*(1-u)*(1-u));
}

double dfuncY(vertex x, double u, double t){
    return dfuncX(x,u,t)*(1-5*(1-u)*(1-u)) -10*(u-1)*funcX(x,u,t);
}

// ======== Diffusion =======================
double diffFunc(double u){
    return u;
}

double dDiffFunc(double u){
    return 1;
}
