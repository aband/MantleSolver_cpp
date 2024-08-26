#include "eutectic.h"

using namespace EUTECTIC;

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

    ste_ = cp_.solid1 * DT_ / L_;

    CD1_ = 0;
    HD1_ = 0;
    CD2s_ = 0;
    HD2s_ = 1;
    CD2l_ = 0;
    HD2l_ = rho_.fluid/rho_.solid1* (1.0/ste_ + cp_.fluid/cp_.solid1);
    CD3s_ = rho_.solid2/rho_.solid1 * Xe_ / 
           (rho_.solid2/rho_.solid1 + Xe_ * (1-rho_.solid2/rho_.solid1)); 
    HD3s_ = 0;
    CD3l_ = rho_.fluid/rho_.solid1 * Xe_;
    HD3l_ = rho_.fluid/rho_.solid1 / ste_;

    nu_ = rho_.solid2/rho_.solid1 + Xe_ * (1-rho_.solid2/ rho_.solid1);
}

// Lines separates different phase regions
double phaseState::HII_III_(const double& CD){

    return 0;
}

double phaseState::HIII_IV_(const double& CD){

    return CD / Xe_/ ste_;
}

double phaseState::HIII_(const double& CD){

    double tmp = (rho_.solid2/rho_.solid1 * rho_.fluid/rho_.solid1) /
                 (ste_ * (rho_.fluid/rho_.solid1 * nu_ - rho_.solid2/rho_.solid1));

    tmp = tmp* (nu_ * CD / (rho_.solid2/rho_.solid1 * Xe_) - 1);

    return tmp;
}

double phaseState::HIV_V_(const double& CD){

    double tmp = cp_.fluid/cp_.solid1 * (1 - CD/(rho_.fluid/rho_.solid1 * Xe_)); 

    tmp = (tmp  + 1.0/ste_) * rho_.fluid/rho_.solid1;

    return tmp;
}

double phaseState::CIII_(const double& HD){

    double tmp = rho_.solid2/rho_.solid1 * Xe_ / nu_ * (1- ste_*HD/ (rho_.fluid/rho_.solid1));

    tmp += ste_*Xe_ * HD;

    return tmp;
}

int phaseState::EvalPhaseRegion(const double& CD,
                                const double& HD){

    // Evaluate current phase state with given information 
    int state = 0; 

    if (CD < std::numeric_limits<double>::epsilon() && HD < HD2s_){
        // Single phase sub-solidus region
        state = 1;    
    } else if (CD > std::numeric_limits<double>::epsilon() && CD < CD3s_ && HD < HII_III_(CD)){
        // Two phase sub-solidus region no liquid
        state = 2;
    } else if ((CD > std::numeric_limits<double>::epsilon() && CD < CD3l_ && 
                HD > std::numeric_limits<double>::epsilon() && HD < HIII_IV_(CD) && HD > 0) || 
               (CD > CD3l_ && CD < CD3s_ && HD < HIII_(CD) && HD > 0)){
        // Eutectic three phases coexists
        state = 3;
    } else if (CD < CD3l_ && HD > HIII_IV_(CD) && HD < HIV_V_(CD)){
        // Super eutectic two-phase region
        state = 4;
    } else if (CD < CD3l_ && HD > HIV_V_(CD)){
        // Single phase super liquidus region
        state = 5;
    }

    return state;
}

int phaseState::EvalPhase(const int& state, 
                          const double& CD,
                          const double& HD){
    double a,b,c;

    switch(state){
        case 1:
            phi.solid1 = 1;
            phi.solid2 = 0;
            phi.fluid  = 0;
            TD         = HD;
        break;

        case 2:
            phi.solid1 = 1-CD/(rho_.solid2/rho_.solid1);
            phi.solid2 = 1-phi.solid1;
            phi.fluid  = 0;
            TD         = HD/ ((1-phi.solid2) + 
                              phi.solid2 * rho_.solid2/rho_.solid1 * cp_.solid2/cp_.solid1);

        break;

        case 3:
            phi.solid2 = (CD - ste_*Xe_*HD) / (rho_.solid2/rho_.solid1);
            phi.fluid  = (ste_*HD/ (rho_.fluid/rho_.solid1)); 
            phi.solid1 = 1-phi.fluid-phi.solid2;

            TD         = 0;

        break;

        case 4:
            a = rho_.fluid/rho_.solid1 * Xe_;
            b = ((1-rho_.fluid/rho_.solid1 * cp_.fluid/cp_.solid1) * CD - 
                        (1- HD) * rho_.fluid/rho_.solid1 * Xe_);
            c = rho_.fluid/rho_.solid1 * (Xe_ * HD - CD/ste_); 

            TD         = (-1*b - sqrt(b*b - 4*a*c)) / (2*a);

            phi.solid2 = 0;
            phi.fluid  = CD/ (rho_.fluid/rho_.solid1 * Xe_ *(1-TD));

            phi.solid1 = 1- phi.solid2 - phi.fluid;
        break;

        case 5:
            phi.solid1 = 0.0;
            phi.solid2 = 0.0;

            phi.fluid  = 1.0;

            TD         = (HD - rho_.fluid/rho_.solid1/ste_) / 
                         (rho_.fluid/rho_.solid1 * cp_.fluid/cp_.solid1);

        break;

        default:

            std::cout << "This is not a valid phase region T_T ." << std::endl;

    }
    return 1;
} 

