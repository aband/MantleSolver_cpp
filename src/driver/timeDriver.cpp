#include "driver.h"
#include "advectiveFlux.h"
#include "diffusiveFlux.h"
#include "edgeFlux.h"
#include "lagrange_tmp.h"
// My RK time stepping solver

int Driver::UpdateFluxAll(const bool& event,
                          const double& dt,
                          Vec * globalhd,
                          Vec * globalcd,
                          Vec * gfluxHD,
                          Vec * gfluxCD){

    // Update global flux vector

    double time = 0.0;

    Vec gHD = *globalhd;
    Vec gCD = *globalcd;
    Vec fHD = *gfluxHD;
    Vec fCD = *gfluxCD;

    //! Get local vectors for solution
    Vec localc;
    Vec localh;

    PetscCall(DMGetLocalVector(dmu, &localc));
    PetscCall(DMGetLocalVector(dmu, &localh));

    PetscCall(DMGlobalToLocalBegin(dmu, gCD, INSERT_VALUES, localc));
    PetscCall(DMGlobalToLocalEnd(dmu, gCD, INSERT_VALUES, localc));
    PetscCall(DMGlobalToLocalBegin(dmu, gHD, INSERT_VALUES, localh));
    PetscCall(DMGlobalToLocalEnd(dmu, gHD, INSERT_VALUES, localh));

    double ** lc;
    double ** lh;

    DMDAVecGetArray(dmu, localc, &lc);
    DMDAVecGetArray(dmu, localh, &lh);

    //! Hook double pointer to meshInfo object for the calculation of smoothness indicator 
	 //! and non linear weights later

    mi.localCD = lc;
    mi.localHD = lh;

    //! Get local vectors for flux
    Vec localfluxcd;
    Vec localfluxhd;

    PetscCall(DMGetLocalVector(dmu, &localfluxcd));
    PetscCall(DMGetLocalVector(dmu, &localfluxhd));

    PetscCall(DMGlobalToLocalBegin(dmu,fCD,INSERT_VALUES,localfluxcd));
    PetscCall(DMGlobalToLocalEnd(dmu,fCD,INSERT_VALUES,localfluxcd));
    PetscCall(DMGlobalToLocalBegin(dmu,fHD,INSERT_VALUES,localfluxhd));
    PetscCall(DMGlobalToLocalEnd(dmu,fHD,INSERT_VALUES,localfluxhd));

    double ** lfcd;
    double ** lfhd;

    PetscCall(DMDAVecGetArray(dmu, localfluxcd, &lfcd));
    PetscCall(DMDAVecGetArray(dmu, localfluxhd, &lfhd));

    UpdateSmoothnessIndicator();
    UpdateNonlinearWgts();

    if (event){

        eventCount++;
        SolveFlow(maxIter, tolUzawa);
        CreateScatterVec();

        // Add up event and output data
        PrintPhaseEvent();
    
        PrintFlowEvent();
        PrintPressureEvent();

        PrintHDEvent();
        PrintCDEvent();
    }

    // self defined time stepping solver 
    vector<double> edgefluxHD;
    vector<double> edgefluxCD;

    UpdateEdgeFluxAll(edgefluxHD, edgefluxCD, advFlux, advFluxBndry, 
                                              advFlux, advFluxBndry);

    for (int j=mi.MPIlocalCellStart[1]; 
             j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1]; j++){
    for (int i=mi.MPIlocalCellStart[0]; 
             i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0]; i++){

        double fluxHD, fluxCD;

        ComputeCellFlux({i,j}, fluxHD, fluxCD, edgefluxHD, edgefluxCD);

        lc[j][i] -= dt*fluxCD;
        lh[j][i] -= dt*fluxHD;

    }}

    // Restore HD and CD
    // (Update ghost region)
    PetscCall(DMDAVecRestoreArray(dmu, localc, &lc));
    PetscCall(DMDAVecRestoreArray(dmu, localh, &lh));
    PetscCall(DMRestoreLocalVector(dmu, &localc));
    PetscCall(DMRestoreLocalVector(dmu, &localh));

    // Restore global flux
    PetscCall(DMDAVecRestoreArray(dmu, localfluxcd, &lfcd));
    PetscCall(DMDAVecRestoreArray(dmu, localfluxhd, &lfhd));
    PetscCall(DMRestoreLocalVector(dmu, &localfluxcd));
    PetscCall(DMRestoreLocalVector(dmu, &localfluxhd));

    return 1;
}

int Driver::RK(){

    double time = 0.0;

    bool event = true;

    while (time < Tmax){

        if ((int)floor(time/dt) % 10 == 1){
            event = true;
        } else {
            event = false;
        }

        //UpdateFluxAll(event);

        time += dt;
    }

    return 1;
}
