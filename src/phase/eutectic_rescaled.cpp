#include "eutectic_rescaled.h"

EUTECTIC::phase::phase(){

    nu    = 6.5*pow(10,6);
    gamma = 1.0/nu;

    Tm0   = 2000;
    Te0   = 1480;
    dT    = Tm0 - Te0;

    L     = 5*pow(10,5);
    cp    = 1200;

    LD    = L/cp/dT;
    TDm0  = 1.0;
    TDe0  = Te0/dT;

    rho   = 3000;
}

double EUTECTIC::phase::GetTDp(const double& TD, 
                               const double& P) const{
    return TD + gamma * P/dT;
}

double EUTECTIC::phase::GetStaticP(const double& zD,
                                   const double& l0) const{
    return rho*g*zD*l0;
}

// Taking dimensionless enthalpy, dimensionless composition and
// Pressure with dimension as input
// Be careful!
int EUTECTIC::phase::evalPhase(const double& HD,
                               const double& CD,
                               const double& P) const{

    

    return 1;
}
