#include "timestepping.h"

/**
 * Explicit Eurler.
 * Using time stepping object provided by Petsc.
 */
PetscErrorCode Explicit(TS ts, PetscReal time, Vec U, Vec F, void* ctx){

    PetscErrorCode    ierr;
    PetscFunctionBeginUser; 

    Ctx * user = (Ctx*) ctx;
    DM dmu = (DM)user->dmu;

    //! Get local vector
    Vec localu;
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, U, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, U, INSERT_VALUES, localu); 

    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    double ** f;
    DMDAVecGetArray(dmu, F, &f);
  
    //! Hook local u to meshinfo
    user->mi->localVals = lu;

    //! Define advection problem first
    user->trPtr->AssignReconstruction(user->trPtr->advection::reconstMethods,
                                      user->trPtr->advection::wenoLevels);

    user->trPtr->SeparateBoundaryLayer(*(user->mi));

    user->trPtr->UpdateNonLinearWgts(*(user->mi),2);

    //user->trPtr->Check(*(user->mi));

    //! Loop through computational domain
    for (int j=user->mi->MPIlocalCellStart[1]; j<user->mi->MPIlocalCellStart[1] + user->mi->MPIlocalCellSize[1]; j++){
    for (int i=user->mi->MPIlocalCellStart[0]; i<user->mi->MPIlocalCellStart[0] + user->mi->MPIlocalCellSize[0]; i++){
        f[j][i] = -1.0*user->trPtr->advFlux(*(user->mi), {i,j}, time);
    }}

    //! Restore array to local vectors.
    DMDAVecRestoreArray(dmu, F, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}
