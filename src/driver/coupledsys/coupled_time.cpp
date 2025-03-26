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
        //advection.computeWgts(ml, mi, h0, allwgtsHD);
        advection.computeWgts(ml, mi, h0, allwgtsHD, location);

        ml.updatesigma(lCD);
        Tensor<weights> allwgtsCD;
        //advection.computeWgts(ml, mi, h0, allwgtsCD);
        advection.computeWgts(ml, mi, h0, allwgtsCD, location);

        // Solve for velocity
        if (t%10 == 0){
            cout << "Darcy-Stokes system solved at : " << t*dt << endl;
            SolveFlow(maxIter, tolUzawa, allwgtsHD, lHD, allwgtsCD, lCD);
            CreateScatterVec();

            printCellAve(mark, &globalHD, mi, "HD");
            printCellAve(mark, &globalCD, mi, "CD");
            PrintFlowEvent(mark);
            PrintPhaseEvent(mark);
            PrintPressureSerialApprox(mark);
            mark ++;
        }

        getflux(allwgtsHD, lHD, allwgtsCD, lCD, lfHD, lfCD);

        PrintEffVel(mark, 2, allwgtsHD, lHD, allwgtsCD, lCD);

        DMDAVecRestoreArray(dmu, fluxHD, &lfHD);
        DMDAVecRestoreArray(dmu, fluxCD, &lfCD);
        DMDAVecRestoreArray(dmu, localHD, &lHD);
        DMDAVecRestoreArray(dmu, localCD, &lCD);
        DMRestoreLocalVector(dmu, &localHD);
        DMRestoreLocalVector(dmu, &localCD);

        VecAXPY(globalHD, -1*dt, fluxHD);
        VecAXPY(globalCD, -1*dt, fluxCD);

   }

    printCellAve(mark, &globalHD, mi, "HD");
    printCellAve(mark, &globalCD, mi, "CD");

    return 1;
}

int Driver::getFluxAll(Vec * fCD, Vec * fHD, Vec * gCD, Vec * gHD, 
                       int t, double dt, int maxIter, double tolUzawa, int interval){

        Vec fluxHD = *fHD;
        Vec fluxCD = *fCD;

        Vec globHD = *gHD;
        Vec globCD = *gCD;

        Vec localHD, localCD; 
        PetscCall(DMGetLocalVector(dmu, &localHD));
        PetscCall(DMGetLocalVector(dmu, &localCD));

        double ** lHD;
        double ** lCD;
        double ** lfHD;
        double ** lfCD;

        PetscCall(DMGlobalToLocalBegin(dmu, globHD, INSERT_VALUES, localHD));
        PetscCall(DMGlobalToLocalEnd(dmu, globHD, INSERT_VALUES, localHD));

        PetscCall(DMDAVecGetArray(dmu, localHD, &lHD););

        PetscCall(DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localCD));
        PetscCall(DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localCD));

        PetscCall(DMDAVecGetArray(dmu, localCD, &lCD););

        PetscCall(DMDAVecGetArray(dmu, fluxHD, &lfHD));
        PetscCall(DMDAVecGetArray(dmu, fluxCD, &lfCD));

        ml.updatesigma(lHD);
        Tensor<weights> allwgtsHD;
        //advection.computeWgts(ml, mi, h0, allwgtsHD);
        advection.computeWgts(ml, mi, h0, allwgtsHD, location);

        ml.updatesigma(lCD);
        Tensor<weights> allwgtsCD;
        //advection.computeWgts(ml, mi, h0, allwgtsCD);
        advection.computeWgts(ml, mi, h0, allwgtsCD, location);

        if (t%interval == 0){
            cout << "Darcy-Stokes system solved at : " << t *dt << endl;
            SolveFlow(maxIter, tolUzawa, allwgtsHD, lHD, allwgtsCD, lCD);
            CreateScatterVec();
        }

        getflux(allwgtsHD, lHD, allwgtsCD, lCD, lfHD, lfCD);

        DMDAVecRestoreArray(dmu, fluxHD, &lfHD);
        DMDAVecRestoreArray(dmu, fluxCD, &lfCD);
        DMDAVecRestoreArray(dmu, localHD, &lHD);
        DMDAVecRestoreArray(dmu, localCD, &lCD);
        DMRestoreLocalVector(dmu, &localHD);
        DMRestoreLocalVector(dmu, &localCD);

    return 1;
}

int Driver::SSP2RK(double dt, double Tmax, int maxIter, double tolUzawa, int interval){

    int mark = 1;

    int Nt = (int)(Tmax/dt);

    Vec gHD_temp, gCD_temp;
    PetscCall(VecDuplicate(globalHD, &gHD_temp));
    PetscCall(VecDuplicate(globalCD, &gCD_temp));

    PetscCall(VecCopy(globalHD, gHD_temp));
    PetscCall(VecCopy(globalCD, gCD_temp));

    for (int t=0; t<Nt; t++) {

        Vec fHD, fCD; 
        PetscCall(VecDuplicate(globalHD, &fHD));
        PetscCall(VecDuplicate(globalCD, &fCD));

        getFluxAll(&fCD, &fHD, &globalCD, &globalHD, 
                   t, dt, maxIter, tolUzawa, interval);

        VecAXPY(gHD_temp, -1*dt, fHD);
        VecAXPY(gCD_temp, -1*dt, fCD);

        Vec fHD2, fCD2; 
        PetscCall(VecDuplicate(globalHD, &fHD2));
        PetscCall(VecDuplicate(globalCD, &fCD2));

        getFluxAll(&fCD2, &fHD2, &gCD_temp, &gHD_temp, 
                   t, dt, maxIter, tolUzawa, interval);

        VecScale(globalHD, 0.5);
        VecScale(globalCD, 0.5);

        VecAXPY(globalHD, 0.5, gHD_temp);
        VecAXPY(globalCD, 0.5, gCD_temp); 

        VecAXPY(globalHD, -0.5*dt, fHD2);
        VecAXPY(globalCD, -0.5*dt, fCD2); 

        if (t%interval == 0){

            printCellAve(mark, &globalHD, mi, "HD");
            printCellAve(mark, &globalCD, mi, "CD");
            PrintFlowEvent(mark);
            PrintPhaseEvent(mark);
            PrintPressureSerialApprox(mark);

            mark ++;
        }

    }

    mark ++;
    printCellAve(mark, &globalHD, mi, "HD");
    printCellAve(mark, &globalCD, mi, "CD");

    return 1;
}
