#include "timestepping.h"

extern "C"{
#include "output.h"
}

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

    user->trPtr->UpdateSmoothnessIndic(*(user->mi));

//    user->trPtr->advection::UpdateNonLinearWgts(*(user->mi));

    user->trPtr->diffusion::UpdateNonLinearWgts(*(user->mi));

//    user->trPtr->advection::UpdateEdgeFlux(*(user->mi));
 
    user->trPtr->diffusion::UpdateEdgeFlux(*(user->mi));

    //! Loop through computational domain
    for (int j=user->mi->MPIlocalCellStart[1]; j<user->mi->MPIlocalCellStart[1] + user->mi->MPIlocalCellSize[1]; j++){
    for (int i=user->mi->MPIlocalCellStart[0]; i<user->mi->MPIlocalCellStart[0] + user->mi->MPIlocalCellSize[0]; i++){
        //f[j][i] = -1.0*user->trPtr->advection::Flux(*(user->mi), {i,j});
        f[j][i] = user->trPtr->diffusion::Flux(*(user->mi), {i,j});
        //cout << "( " << i << ", " << j << " )" << " Flux : " << f[j][i] << ";  ";
    }}//cout << endl;}

    //! Compute diffusion flux second ======================================
//    user->trPtr->AssignReconstruction(user->trPtr->diffusion::reconstMethodsVert,
//                                      user->trPtr->diffusion::wenoLevelsVert);

 //   user->trPtr->UpdateNonLinearWgts(*(user->mi),2);


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
//    user->trPtr->AssignReconstruction(user->trPtr->advection::reconstMethods,
//                                      user->trPtr->advection::wenoLevels);

//    user->trPtr->AssignBoundaryMethods(user->trPtr->advection::boundaryCells,
//                                       user->trPtr->advection::interiorCells,
//                                       user->trPtr->advection::boundaryLevels,
//                                       user->trPtr->advection::interiorLevels);

    //user->trPtr->UpdateNonLinearWgts(*(user->mi),2);

    //! Get MPI local part of Jacobian matrix
    int rstart, rend;
    MatGetOwnershipRange(J, &rstart, &rend);

    for (int row = rstart; row<rend; row++){
        indice global = Bend(*(user->mi), row);

        unordered_map<int,double> deriv = user->trPtr->advection::derivFlux(*(user->mi), global);

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

    // Output Jacobian
    //std::string tmp = std::to_string(time);

    //DrawMat(J,tmp.c_str());

    //! Restore array to local vectors
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}
