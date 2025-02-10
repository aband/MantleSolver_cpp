#include "advectiveFlux.h"
#include "diffusiveFlux.h"
#include "edgeFlux.h"
#include "driver.h"
#include "lagrange_tmp.h"

PetscErrorCode Explicit(TS ts, PetscReal time, Vec U, Vec F, void* ctx){

    PetscFunctionBeginUser;

    ctx_driver * user = (ctx_driver*) ctx;
    DM dmu = (DM)user->driver->dmu;

    //! Get individual nested vectors 
    Vec C, H;

    PetscCall(VecNestGetSubVec(U, 0, &C));
    PetscCall(VecNestGetSubVec(U, 1, &H));

    //! Get local vectors for solution
    Vec localc;
    Vec localh;

    DMGetLocalVector(dmu, &localc);
    DMGetLocalVector(dmu, &localh);

    //PetscCall(DMGetLocalVector(dmu, C, INSERT_VALUES, localc));
    //PetscCall(DMGetLocalVector(dmu, H, INSERT_VALUES, localh));

    PetscCall(DMGlobalToLocalBegin(dmu, C, INSERT_VALUES, localc));
    PetscCall(DMGlobalToLocalEnd(dmu, C, INSERT_VALUES, localc));
    PetscCall(DMGlobalToLocalBegin(dmu, H, INSERT_VALUES, localh));
    PetscCall(DMGlobalToLocalEnd(dmu, H, INSERT_VALUES, localh));

    double ** lc;
    double ** lh;

    DMDAVecGetArray(dmu, localc, &lc);
    DMDAVecGetArray(dmu, localh, &lh);

    user->driver->mi.localCD = lc;
    user->driver->mi.localHD = lh;

    //! Get local vectors for flux vector
    PetscCall(VecDuplicate(U, &F));

    Vec CF;
    Vec HF;

    PetscCall(VecNestGetSubVec(F, 0, &CF));
    PetscCall(VecNestGetSubVec(F, 1, &HF));

    //! Flux vector does not need to assign ghost region here
    double ** cf;
    double ** hf;

    PetscCall(DMDAVecGetArray(dmu, CF, &cf));
    PetscCall(DMDAVecGetArray(dmu, HF, &hf));

    //! Update smoothness indicator
  //  user->driver->UpdateSmoothnessIndicator();

    //! Update nonlinear weights using new smoothness indicator
  //  user->driver->UpdateNonlinearWgts();

    if ((int)floor(time/user->dt) % 10 == 1){
        // solve for stokes equation when real time meets some standards
        cout << "Darcy-Stokes coupled system solved at time: " << time << endl;

        // Print selected variables
        // Increase internal event count by 1
        user->driver->SolveFlow(user->maxIter, user->tolUzawa);
        user->driver->CreateScatterVec();

        // eventCount++ should always be called before calling print event functions
        //user->driver->eventCount++;
        //user->driver->PrintPhaseEvent();
    
        //user->driver->PrintFlowEvent();
        //user->driver->PrintPressureEvent();

        //user->driver->PrintHDEvent();
        //user->driver->PrintCDEvent();
    }

    //! Compute updated edge flux
    vector<double> edgefluxHD;
    vector<double> edgefluxCD;

  //  user->driver->UpdateEdgeFluxAll(edgefluxHD, edgefluxCD, advFlux, advFluxBndry, advFlux, advFluxBndry);

    for (int j=user->driver->mi.MPIlocalCellStart[1]; 
             j<user->driver->mi.MPIlocalCellStart[1] + user->driver->mi.MPIlocalCellSize[1]; j++){
    for (int i=user->driver->mi.MPIlocalCellStart[0]; 
             i<user->driver->mi.MPIlocalCellStart[0] + user->driver->mi.MPIlocalCellSize[0]; i++){

        double fluxHD, fluxCD;

   //     user->driver->ComputeCellFlux({i,j}, fluxHD, fluxCD, edgefluxHD, edgefluxCD);

        //cf[j][i] = -1.0 * fluxCD;

        //hf[j][i] = -1.0 * fluxHD;

        cf[j][i] = 0.0;

        hf[j][i] = 100.0;

    }}

    // Restore flux vector
    PetscCall(DMDAVecRestoreArray(dmu, CF, &cf));
    PetscCall(DMDAVecRestoreArray(dmu, HF, &hf));

    // Restore flux 
    PetscCall(DMDAVecRestoreArray(dmu, localc, &lc));
    PetscCall(DMDAVecRestoreArray(dmu, localh, &lh));
    PetscCall(DMRestoreLocalVector(dmu, &localc));
    PetscCall(DMRestoreLocalVector(dmu, &localh));

    VecView(F, PETSC_VIEWER_STDOUT_WORLD);
    VecView(U, PETSC_VIEWER_STDOUT_WORLD);

    PetscFunctionReturn(0);
}
