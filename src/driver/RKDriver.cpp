#include "driver.h"
// My RK time stepping solver

int Driver::RK(const double& dt, 
               const double& Tmax,
               const double& tolUzawa,
               const int& maxIter){

    double time = 0.0;

    //! Get local vectors for solution
    Vec localc;
    Vec localh;

    DMGetLocalVector(dmu, &localc);
    DMGetLocalVector(dmu, &localh);

    PetscCall(DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localc));
    PetscCall(DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localc));
    PetscCall(DMGlobalToLocalBegin(dmu, globalHD, INSERT_VALUES, localh));
    PetscCall(DMGlobalToLocalEnd(dmu, globalHD, INSERT_VALUES, localh));

    double ** lc;
    double ** lh;

    DMDAVecGetArray(dmu, localc, &lc);
    DMDAVecGetArray(dmu, localh, &lh);

    mi.localCD = lc;
    mi.localHD = lh;

    // self defined time stepping solver 
    while(time < Tmax) {


        time += dt;
    }

    return 1;
}
