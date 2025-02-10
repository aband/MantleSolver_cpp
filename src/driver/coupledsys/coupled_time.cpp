#include "driver.h"
#include "print.h"

int Driver::RK(double dt, double Tmax, int maxIter, double tolUzawa){

    int mark = 1;

    int Nt = (int)(Tmax/dt);

    printCellAve(mark, &globalHD, mi, "HD");
    printCellAve(mark, &globalCD, mi, "CD");

    for (int t=0; t<Nt; t++) {

        Vec localHD, localCD; 
        PetscCall(DMGetLocalVector(dmu, &localHD));
        PetscCall(DMGetLocalVector(dmu, &localCD));
    
        Vec fluxHD, fluxCD;
        PetscCall(VecDuplicate(globalHD, &fluxHD));
        PetscCall(VecDuplicate(globalCD, &fluxCD));

        double ** lHD;
        double ** lCD;
        double ** lfHD;
        double ** lfCD;

        PetscCall(DMGlobalToLocalBegin(dmu, globalHD, INSERT_VALUES, localHD));
        PetscCall(DMGlobalToLocalEnd(dmu, globalHD, INSERT_VALUES, localHD));

        PetscCall(DMDAVecGetArray(dmu, localHD, &lHD););

        PetscCall(DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localCD));
        PetscCall(DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localCD));

        PetscCall(DMDAVecGetArray(dmu, localCD, &lCD););

        PetscCall(DMDAVecGetArray(dmu, fluxHD, &lfHD));
        PetscCall(DMDAVecGetArray(dmu, fluxCD, &lfCD));


        ml.updatesigma(lHD);
        Tensor<weights> allwgtsHD;
        advection.computeWgts(ml, mi, h0, allwgtsHD);

        ml.updatesigma(lCD);
        Tensor<weights> allwgtsCD;
        advection.computeWgts(ml, mi, h0, allwgtsCD);

        // Solve for velocity
        //if (t  == 0){
            cout << "Darcy-Stokes system solved at : " << t*dt << endl;
            SolveFlow(maxIter, tolUzawa, allwgtsHD, lHD, allwgtsCD, lCD);
            CreateScatterVec();
           //PrintFlowEvent(mark);
           //PrintPhaseEvent(mark);
        //}
/*
        if (t%2 == 0){
             PrintEffVel(mark-1, 2, allwgtsHD, lHD, allwgtsCD, lCD);
        }
*/
        getflux(allwgtsHD, lHD, allwgtsCD, lCD, lfHD, lfCD);

        DMDAVecRestoreArray(dmu, fluxHD, &lfHD);
        DMDAVecRestoreArray(dmu, fluxCD, &lfCD);
        DMDAVecRestoreArray(dmu, localHD, &lHD);
        DMDAVecRestoreArray(dmu, localCD, &lCD);
        DMRestoreLocalVector(dmu, &localHD);
        DMRestoreLocalVector(dmu, &localCD);

        VecAXPY(globalHD, -1*dt, fluxHD);
        VecAXPY(globalCD, -1*dt, fluxCD);

        if (t%5 == 0){
           printCellAve(mark, &globalHD, mi, "HD");
           printCellAve(mark, &globalCD, mi, "CD");
           PrintFlowEvent(mark);
           PrintPhaseEvent(mark);

           mark ++;
       }
    }

    printCellAve(mark, &globalHD, mi, "HD");
    printCellAve(mark, &globalCD, mi, "CD");

    return 1;
}
