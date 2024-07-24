#include "simplePhase.h"

phaseState::phaseState(){
    // A construtor initialize constant physical attributes

    L_  = 4e5;  // Latent heat in melting 
    Xe_ = 0.7;  // Eutectic composition of solid 2
    T1_ = 1350; // Melting temperature for solid 1
                // In Eutectic phase behavior, melting temperature 1
                // belongs to solid 1 that has lower melting temperature
    Te_ = 1227; // Eutectic temperature

    // density
    rho_.solid1 = 3e3;
    rho_.solid2 = 3e3;
    rho_.fluid  = 3e3;

    // specific enthalpy
    cp_.solid1  = 1200;
    cp_.solid2  = 1200;
    cp_.fluid   = 1200;

    // thermal conductivity
    kappa_.solid1 = 5.2;
    kappa_.solid2 = 4.7;
    kappa_.fluid  = 4.7;

    // Reference physical attributes for nondimensionalization
    DT_ = T1_ - Te_;
    hc_ = cp_.solid1 * DT_;
    Hc_ = hc_ * rho_.solid1;
    Cc_ = rho_.solid1;

    ste_ = cp_.solid1 * DT_ / L;

    CD1_ = 0;
    HD1_ = 0;
    CD2s_ = 0;
    HD2s_ = 1;
    CD2l_ = 0;
    HD2l_ = rho_.fluid/rho_.solid1* (1.0/ste_ + cp_.fluid/cp_.solid1);
    CD3s_ = rho_.solid2/rho_.sold1 * Xe_ / 
           (rho_.solid2/rho.solid1 + Xe_ * (1-rho_.solid2/rho.solid1)); 
    HD3s_ = 0;
    CD3l_ = rho_.fluid/rho_.solid1 * Xe_;
    HD3l_ = rho_.fluid/rho_.solid1 / ste_;
}

double 

int phaseState::EvalPhase(double CD, double HD){

    // Evaluate current phase state with given information 
    int state = 0; 

    if (CD < 2*DBL_EPSILON && HD < HD2s_){
        // Single phase sub-solidus region
        state = 1;    
    } else if (CD > 2*DBL_EPSILON && CD < CD3s_ && HD < 0){
        // Two phase sub-solidus region no liquid
        state = 2;
    } else if ((CD > 2*DBL_EPSILON && CD < CD3l_ && 
                HD > 2*DBL_EPSILON && HD < CD/Xe_/ste_) || 
               (CD > CD3l_ && CD < CD3s_)){


    }


    return state;
}


inline double getSysD(Phase * phase, Phi * phi){

    double sysD = phi->ice * phase->kappa_ii+ 
                  phi->sal * phase->kappa_si+ 
                  phi->bri * phase->kappa_bi;

    return sysD;
}

int phaseState::EvalVolumeFrac(int state){

    switch(){
        case 1:

        break;
    }

} 

