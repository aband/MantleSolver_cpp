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

    //! Compute advection flux first =======================================
    user->trPtr->AssignReconstruction(user->trPtr->advection::reconstMethods,
                                      user->trPtr->advection::wenoLevels);

    user->trPtr->AssignBoundaryMethods(user->trPtr->advection::boundaryCells,
                                       user->trPtr->advection::interiorCells,
                                       user->trPtr->advection::boundaryLevels,
                                       user->trPtr->advection::interiorLevels);

    user->trPtr->UpdateNonLinearWgts(*(user->mi),2);

    //user->trPtr->Check(*(user->mi));

    //! Loop through computational domain
    for (int j=user->mi->MPIlocalCellStart[1]; j<user->mi->MPIlocalCellStart[1] + user->mi->MPIlocalCellSize[1]; j++){
    for (int i=user->mi->MPIlocalCellStart[0]; i<user->mi->MPIlocalCellStart[0] + user->mi->MPIlocalCellSize[0]; i++){
        f[j][i] = -1.0*user->trPtr->advFlux(*(user->mi), {i,j}, time);
    }}

    //! Compute diffusion flux second ======================================



    //! Restore array to local vectors.
    DMDAVecRestoreArray(dmu, F, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}

/**
 * Compute jacobian requiared for implicit time stepping
 */
PetscErrorCode FormJacobian(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void* ctx){
    PetscErrorCode    ierr;
    PetscFunctionBeginUser;

    Ctx * user = (Ctx*)ctx;
    DM dmu = (DM)user->dmu;

    //! Get local vector
    Vec localu;
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, U, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, U, INSERT_VALUES, localu); 

    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    //! Define advection problem first
    user->trPtr->AssignReconstruction(user->trPtr->advection::reconstMethods,
                                      user->trPtr->advection::wenoLevels);

    user->trPtr->AssignBoundaryMethods(user->trPtr->advection::boundaryCells,
                                       user->trPtr->advection::interiorCells,
                                       user->trPtr->advection::boundaryLevels,
                                       user->trPtr->advection::interiorLevels);

    //user->trPtr->UpdateNonLinearWgts(*(user->mi),2);

    //! Get MPI local part of Jacobian matrix
    int rstart, rend;
    MatGetOwnershipRange(J, &rstart, &rend);

    for (int row = rstart; row<rend; row++){
        indice global = Bend(*(user->mi), row);

        unordered_map<int,double> deriv = user->trPtr->derivAdvFlux(*(user->mi), global, time);

        for (auto & dVal: deriv){
            ierr = MatSetValue(J,row,dVal.first,-1.0*dVal.second,INSERT_VALUES);CHKERRQ(ierr);
        }

    }

    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (J != Jp){
        ierr = MatAssemblyBegin(Jp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Jp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

    //MatView(J,PETSC_VIEWER_STDOUT_WORLD);

    //! Restore array to local vectors
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}
