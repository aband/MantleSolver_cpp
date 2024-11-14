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

    PetscCall(DMGetLocalVector(dmu, C, INSERT_VALUES, localc));
    PetscCall(DMGetLocalVector(dmu, H, INSERT_VALUES, localh));

    double ** lc;
    double ** lh;

    DMDAVecGetArray(dmu, localc, &lc);
    DMDAVecGetArray(dmu, localh, &lh);

    user->mi->localCD = lc;
    user->mi->localHD = lh;

    //! Get local vectors for flux vector
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
    user->mlpPtr->UpdateSmoothnessIndic(user->mi, user->mi->localCD, "CD");
    user->mlpPtr->UpdateSmoothnessIndic(user->mi, user->mi->localHD, "HD");

    //! Update nonlinear weights using new smoothness indicator
    user->mluseAdv->UpdateNonLinearWgts(mi,user->locpack->locSet,
                                           user->locpack->funcSet,
                                           user->locpack->fieldNames);

    if (){
        // solve for stokes equation when real time meets some standards

    }

    for (int j=user->mi->MPIlocalCellStart[1]; j<user->mi->MPIlocalCellStart[1] + user->mi->MPIlocalCellSize[1]; j++){
    for (int i=user->mi->MPIlocalCellStart[0]; i<user->mi->MPIlocalCellStart[0] + user->mi->MPIlocalCellSize[0]; i++){
 
        cf[j][i] = -1.0 * user->

    }}

    // Restore flux vector
    PetscCall(DMDAVecRestoreArray(dmu, CF, &cf));
    PetscCall(DMDAVecRestoreArray(dmu, HF, &hf));

    // Restore flux 
    PetscCall(DMDAVecRestoreArray(dmu, localc, &lc));
    PetscCall(DMRestoareLocalVector(dmu, &localx));

    PetscFunctionReturn(0);
}
