#include "myFunc.h"
#include "param.h"
#include <petsc.h>

// Assign physical properties
void AssignPhyProperties(PhysProperty * pp){

    pp->theta = 0.0;     //Additional parameter accounting for singularity
    pp->mu_s  = 1e19;    //Viscosity solid
    pp->mu_f  = 1.0;     //Viscosity fluid
    pp->rho_f = 2800;    //Density fluid (should be defined later)
    pp->rho_s = 3300;    //Density solid (should be defined later)
    pp->gx    = 0.0;     //Gravitational force in x direction (of course 0)
    pp->gy    = -10.0;   //Gravitational force in y direction ();
    pp->invk0 = 1.0/(1e-8); //Inverse of permeability
    pp->phi0  = 0.4;        //Initial homogenous porosity (if constant)
    pp->U0    = 1e-9;       //Specific velocity
    pp->L0    = 160*1000;   //Specific length (physical domain size)

    // Density of two species in the mixture
    pp->rho_1 = 3.25*1000; // kg/m^3 (olivine)
    pp->rho_2 = 1600;      // kg/m^3 (rock)

    // Non dimensionalization parameters

    double rho_r = pp->rho_f*pp->phi0 + 
                   pp->rho_s*(1-pp->phi0); 

    pp->l0    = pow(pp->mu_s/pp->invk0/pp->mu_f,0.5); //Nondim length parameter
    pp->p0    = pp->gy*pp->l0*rho_r;                  //Nondim pressure parameter
    pp->u0    = pp->gy*rho_r/pp->mu_f/pp->invk0;      //Nondim velocity parameter

    pp->l = 20/pp->l0;
}

// Assign initial porosity
// Only used initially
double AssignPorosity(const vertex& point, PhysProperty * pp){
    if (abs(point[1]) < 120*1000/pp->l0 && abs(point[0]) < abs(point[1]) + pp->l){
        double value = 0.05*pow((120*1000/pp->l0 - abs(point[1]))/(120*1000/pp->l0),2) * 
                               (1-abs(point[0])/(abs(point[1])+pp->l));

        return value;
    } else {
        return 0.0;
    }
}

// Porosity of solid
double AssignPorosity(double phi_f){

    return 1-phi_f;
}

// ================================================================================

inline double InitCD(const vertex& point, PhysProperty * pp){

    if (abs(point[1]) < 120*1000/pp->l0 && abs(point[0]) < abs(point[1]) + pp->l){

        return 0.4;
    } else {
        return 0.2;
    }

}

inline double InitHD(const vertex& point, PhysProperty * pp){

    if (abs(point[1]) < 120*1000/pp->l0 && abs(point[0]) < abs(point[1]) + pp->l){

        // Above Eutectic
        return 1.0;
    } else {
        // Below eutectic
        return -0.05;
    }

}

double ComputePorosity(const vertex& point, Phase * phase){

    double CD = InitCD(point, phase->pp);
    double HD = InitHD(point, phase->pp);

    int state = phase->pPtr->EvalPhaseRegion(CD,HD);
    phase->pPtr->EvalPhase(state, CD, HD);

    return phase->pPtr->phi.fluid;
}

int PorosityOut(double xstart, double ystart, double L, double H, int seed,
                Phase * phase){

    FILE * fp = fopen("InitPoro.dat","w");

    double hx = L/(double)seed;
    double hy = H/(double)seed;

    //phaseState * pPtr = new phaseState();

    for (int j=0; j<seed; j++){
    for (int i=0; i<seed; i++){
        vertex point {xstart + hx*i , ystart + hy*j}; 
        fprintf(fp, "%f ", ComputePorosity(point,phase));
    }fprintf(fp, "\n");}

    fclose(fp);

    return 1;

}
