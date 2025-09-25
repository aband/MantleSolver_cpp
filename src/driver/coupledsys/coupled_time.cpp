#include "driver.h"
#include "print.h"

int Driver::RK(double dt, double Tmax, int maxIter, double tolUzawa){

    int mark = 1 + start;

    int Nt = (int)(Tmax/dt);

    //printCellAve(mark, &globalHD, mi, "HD");
    //printCellAve(mark, &globalCD, mi, "CD");

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
        //cout << t*dt << endl;
        //if (t == 0){
            cout << "Darcy-Stokes system solved at : " << t*dt << endl;
            SolveFlow(maxIter, tolUzawa, allwgtsHD, lHD, allwgtsCD, lCD);
            CreateScatterVec();
       //}

        getflux(allwgtsHD, lHD, allwgtsCD, lCD, lfHD, lfCD);

        if (t % 200 ==0){
        PrintEffVel(mark, 2, allwgtsHD, lHD, allwgtsCD, lCD);
        }

        DMDAVecRestoreArray(dmu, fluxHD, &lfHD);
        DMDAVecRestoreArray(dmu, fluxCD, &lfCD);
        DMDAVecRestoreArray(dmu, localHD, &lHD);
        DMDAVecRestoreArray(dmu, localCD, &lCD);
        DMRestoreLocalVector(dmu, &localHD);
        DMRestoreLocalVector(dmu, &localCD);

        if (t%200 == 0){
        printCellAve(mark, &globalHD, mi, "HD");
        printCellAve(mark, &globalCD, mi, "CD");
        //PrintFlowEvent(mark);
        PrintPhaseEvent(mark);
        PrintPressureSerialApprox(mark);

        mark ++;
		  cout << "Output mark = " << mark << endl;
        }

        VecAXPY(globalHD, -1*dt, fluxHD);
        VecAXPY(globalCD, -1*dt, fluxCD);
   }

    //printCellAve(mark, &globalHD, mi, "HD");
    //printCellAve(mark, &globalCD, mi, "CD");

    return 1;
}

int Driver::velocitycamera(Tensor<vertexSet>& phasevel_vert, 
                           Tensor<vertexSet>& phasevel_hori, 
                           Tensor<vertexSet>& effvel_vert, 
                           Tensor<vertexSet>& effvel_hori, 
                           Tensor<vertexSet>& solidvel_vert, 
                           Tensor<vertexSet>& solidvel_hori, 
                           Vec * gCD, Vec * gHD, 
                           int t, double dt, int maxIter, 
                           double tolUzawa){

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

        ml.updatesigma(lHD);
        Tensor<weights> allwgtsHD;
        //advection.computeWgts(ml, mi, h0, allwgtsHD);
        advection.computeWgts(ml, mi, h0, allwgtsHD, location);

        ml.updatesigma(lCD);
        Tensor<weights> allwgtsCD;
        //advection.computeWgts(ml, mi, h0, allwgtsCD);
        advection.computeWgts(ml, mi, h0, allwgtsCD, location);

        cout << "Darcy-Stokes system solved at : " << t *dt << endl;
        SolveFlow(maxIter, tolUzawa, allwgtsHD, lHD, allwgtsCD, lCD);
        CreateScatterVec();

        updateVel_Pause(phasevel_vert, phasevel_hori, 
                        effvel_vert,   effvel_hori,
                        solidvel_vert, solidvel_hori,
                        allwgtsHD, lHD, allwgtsCD, lCD);

        //PrintEffVel(1, 2, allwgtsHD, lHD, allwgtsCD, lCD);

        DMDAVecRestoreArray(dmu, localHD, &lHD);
        DMDAVecRestoreArray(dmu, localCD, &lCD);
        DMRestoreLocalVector(dmu, &localHD);
        DMRestoreLocalVector(dmu, &localCD);

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

int Driver::getFluxAll(const Tensor<vertexSet>& phasevel_vert, 
                       const Tensor<vertexSet>& phasevel_hori, 
                       const Tensor<vertexSet>& effvel_vert, 
                       const Tensor<vertexSet>& effvel_hori, 
                       const Tensor<vertexSet>& solidvel_vert, 
                       const Tensor<vertexSet>& solidvel_hori,
                       Vec * fCD, Vec * fHD, Vec * gCD, Vec * gHD, 
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

        getflux(phasevel_vert, phasevel_hori,
                effvel_vert,   effvel_hori,
                solidvel_vert, solidvel_hori,
                allwgtsHD, lHD, allwgtsCD, lCD, lfHD, lfCD); 

        DMDAVecRestoreArray(dmu, fluxHD, &lfHD);
        DMDAVecRestoreArray(dmu, fluxCD, &lfCD);
        DMDAVecRestoreArray(dmu, localHD, &lHD);
        DMDAVecRestoreArray(dmu, localCD, &lCD);
        DMRestoreLocalVector(dmu, &localHD);
        DMRestoreLocalVector(dmu, &localCD);

    return 1;
}

int Driver::SSP2RK(double dt, double Tmax, int maxIter, double tolUzawa, int interval){

    int mark = 1 + start;

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

        if (t%1 == 0){

            printCellAve(mark, &globalHD, mi, "HD");
            printCellAve(mark, &globalCD, mi, "CD");
            //PrintFlowEvent(mark);
            PrintPhaseEvent(mark);
            PrintPressureSerialApprox(mark);

            mark ++;
        }

    }

    //mark ++;
    //printCellAve(mark, &globalHD, mi, "HD");
    //printCellAve(mark, &globalCD, mi, "CD");

    return 1;
}

int Driver::RK_Pause(double dt, double Tmax, int maxIter, 
                     double tolUzawa, int interval){

    int mark = 1 + start;

    int Nt = (int)(Tmax/dt);

    Tensor<vertexSet> effvel_vert = Tensor<vertexSet>(2);
    effvel_vert.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});
    Tensor<vertexSet> effvel_hori = Tensor<vertexSet>(2);
    effvel_hori.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<vertexSet> phasevel_vert = Tensor<vertexSet>(2);
    phasevel_vert.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});
    Tensor<vertexSet> phasevel_hori = Tensor<vertexSet>(2);
    phasevel_hori.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<vertexSet> solidvel_vert = Tensor<vertexSet>(2);
    solidvel_vert.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});
    Tensor<vertexSet> solidvel_hori = Tensor<vertexSet>(2);
    solidvel_hori.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    for (int t=0; t<Nt; t++){

        Vec fHD, fCD; 
        PetscCall(VecDuplicate(globalHD, &fHD));
        PetscCall(VecDuplicate(globalCD, &fCD));
 
        if (t% interval == 0){
            velocitycamera(phasevel_vert, phasevel_hori, 
                           effvel_vert,   effvel_hori,
                           solidvel_vert, solidvel_hori,
                           &globalCD, &globalHD, 
                           t, dt, maxIter, tolUzawa);

            PrintTensorVel(phasevel_hori, mark, "phasevel_pause");
            PrintTensorVel(effvel_hori, mark, "effvel_pause");
            PrintTensorVel(solidvel_hori, mark, "solidvel_pause");
        }

        getFluxAll(phasevel_vert, phasevel_hori, 
                   effvel_vert,   effvel_hori,
                   solidvel_vert, solidvel_hori,
                   &fCD, &fHD, &globalCD, &globalHD, 
                   t, dt, maxIter, tolUzawa, interval);

        VecAXPY(globalHD, -1*dt, fHD);
        VecAXPY(globalCD, -1*dt, fCD);

        printCellAve(mark, &globalHD, mi, "HD");
        printCellAve(mark, &globalCD, mi, "CD");
        PrintPhaseEvent(mark);
        mark++;
    }

    mark++;
    printCellAve(mark, &globalHD, mi, "HD");
    printCellAve(mark, &globalCD, mi, "CD");

    return 1;
}

int Driver::SSP2RK_Pause(double dt, double Tmax, int maxIter, 
                         double tolUzawa, int interval){

    int mark = 1 + start;

    int Nt = (int)(Tmax/dt);

    Tensor<vertexSet> effvel_vert = Tensor<vertexSet>(2);
    effvel_vert.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});
    Tensor<vertexSet> effvel_hori = Tensor<vertexSet>(2);
    effvel_hori.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<vertexSet> phasevel_vert = Tensor<vertexSet>(2);
    phasevel_vert.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});
    Tensor<vertexSet> phasevel_hori = Tensor<vertexSet>(2);
    phasevel_hori.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<vertexSet> solidvel_vert = Tensor<vertexSet>(2);
    solidvel_vert.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});
    Tensor<vertexSet> solidvel_hori = Tensor<vertexSet>(2);
    solidvel_hori.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Vec gHD_temp, gCD_temp;
    PetscCall(VecDuplicate(globalHD, &gHD_temp));
    PetscCall(VecDuplicate(globalCD, &gCD_temp));

    PetscCall(VecCopy(globalHD, gHD_temp));
    PetscCall(VecCopy(globalCD, gCD_temp));

    int mark2 = 1 + start;
    for (int t=0; t<Nt; t++){

        Vec fHD, fCD; 
        PetscCall(VecDuplicate(globalHD, &fHD));
        PetscCall(VecDuplicate(globalCD, &fCD));

        if (t% interval == 0){
            velocitycamera(phasevel_vert, phasevel_hori, 
                           effvel_vert,   effvel_hori,
                           solidvel_vert, solidvel_hori,
                           &globalCD, &globalHD, 
                           t, dt, maxIter, tolUzawa);

            PrintFlowEvent(mark2);
            PrintTensorVel(phasevel_hori, mark2, "phasevel_pause");
            PrintTensorVel(effvel_hori, mark2, "effvel_pause");
            PrintTensorVel(solidvel_hori, mark2, "solidvel_pause");
            mark2++;
        }

        getFluxAll(phasevel_vert, phasevel_hori, 
                   effvel_vert,   effvel_hori,
                   solidvel_vert, solidvel_hori,
                   &fCD, &fHD, &globalCD, &globalHD, 
                   t, dt, maxIter, tolUzawa, interval);

        VecAXPY(gHD_temp, -1*dt, fHD);
        VecAXPY(gCD_temp, -1*dt, fCD);

        Vec fHD2, fCD2; 
        PetscCall(VecDuplicate(globalHD, &fHD2));
        PetscCall(VecDuplicate(globalCD, &fCD2));

        getFluxAll(phasevel_vert, phasevel_hori, 
                   effvel_vert,   effvel_hori,
                   solidvel_vert, solidvel_hori,
                   &fCD2, &fHD2, &gCD_temp, &gHD_temp, 
                   t, dt, maxIter, tolUzawa, interval);

        VecScale(globalHD, 0.5);
        VecScale(globalCD, 0.5);

        VecAXPY(globalHD, 0.5, gHD_temp);
        VecAXPY(globalCD, 0.5, gCD_temp); 

        VecAXPY(globalHD, -0.5*dt, fHD2);
        VecAXPY(globalCD, -0.5*dt, fCD2); 

        printCellAve(mark, &globalHD, mi, "HD");
        printCellAve(mark, &globalCD, mi, "CD");
        PrintPhaseEvent(mark);
        mark++;
    }

    return 1;
}
