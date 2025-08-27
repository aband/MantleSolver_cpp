#include "driver.h"

int Driver::rk1(){

    int mark = 1+start;

    int Nt = (int)(Tmax/dt);

    for (int t=0; t<Nt; t++){

        Vec fluxHD;
        Vec fluxCD;

        VecDuplicate(globalHD, &fluxHD);
        VecDuplicate(globalCD, &fluxCD);

        getfluxall(&fluxHD, &fluxCD, true, dt*Nt);

        VecAXPY(globalHD, -1*dt, fluxHD);
        VecAXPY(globalCD, -1*dt, fluxCD);
    }

    return 1;
}
