#include "advectiveFlux.h"
#include "diffusiveFlux.h"
#include ""

PetscErrorCode Explicit(TS ts, PetscReal time, Vec U, Vec F, void* ctx){

    PetscFunctionBeginUser;

    Ctx * user = (Ctx*) ctx;
    DM dmu = (DM)user->dmu;

    //! Get individual nested vectors 
    Vec C, H;

    PetscCall(VecNestGetSubVec(U, 0, &C));
    PetscCall(VecNestGetSubVec(U, 1, &H));

    //! Get local vectors
    Vec localc;
    Vec localh;

    DMGetLocalVector(dmu, C, INSERT_VALUES, localc);
    DMGetLocalVector(dmu, H, INSERT_VALUES, localh);

    double ** lc;
    double ** lh;

    DMDAVecGetArray(dmu, localc, &lc);
    DMDAVecGetArray(dmu, localh, &lh);

    user->mi->local


    PetscFunctionReturn(0);
}
