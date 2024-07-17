#include "simplePhase.h"

inline double getSysD(Phase * phase, Phi * phi){

    double sysD = phi->ice * phase->kappa_ii+ 
                  phi->sal * phase->kappa_si+ 
                  phi->bri * phase->kappa_bi;

    return sysD;
}



int SetEutecticPhase(Phase * phase){

    phase->L  = 4e5;
    phase->Xe = 0.7;
    phase->T1 = 1350;
    phase->Te = 1227;

    phase->DT = phase->T1-phase->Te;

    phase->rho_ice = 3e3;
    phase->rho_sal = 3e3;
    phase->rho_bri = 3e3;

    phase->rho_bi  = phase->rho_bri/ phase->rho_ice;
    phase->rho_si  = phase->rho_sal/ phase->rho_ice;

    phase->cp_ice = 1200;
    phase->cp_sal = 1200;
    phase->cp_bri = 1200;

    phase->cp_bi  = phase->cp_bri/ phase->cp_ice;
    phase->cp_si  = phase->cp_sal/ phase->cp_ice;

    phase->Ste = phase->cp_ice*phase->DT/phase->L; 

    // Constant thermal conductivity right now
    phase->kappa_ice = 5.2;
    phase->kappa_sal = 4.7;
    phase->kappa_bri = 4.7;

    phase->kappa_ii = phase->kappa_ice/phase->kappa_ice;
    phase->kappa_si = phase->kappa_sal/phase->kappa_ice;
    phase->kappa_bi = phase->kappa_bri/phase->kappa_ice;

    return 0;
}
