#include "timestepping.h"

extern "C"{
#include "output.h"
}

/**
 * Explicit advection-diffusion.
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

    // Update smoothness indicators for all reconstruction levels altogether
    user->trPtr->UpdateSmoothnessIndic(*(user->mi));

    // Update non linear weights based on updated smoothness indicators
    user->trPtr->advection::UpdateNonLinearWgts(*(user->mi));
    user->trPtr->diffusion::UpdateNonLinearWgts(*(user->mi));

    // Update flux defined on each edge
    user->trPtr->advection::UpdateEdgeFlux(*(user->mi));
    user->trPtr->diffusion::UpdateEdgeFlux(*(user->mi));

    //! Loop through computational domain
    for (int j=user->mi->MPIlocalCellStart[1]; j<user->mi->MPIlocalCellStart[1] + user->mi->MPIlocalCellSize[1]; j++){
    for (int i=user->mi->MPIlocalCellStart[0]; i<user->mi->MPIlocalCellStart[0] + user->mi->MPIlocalCellSize[0]; i++){
        f[j][i] = -1.0*user->trPtr->advection::Flux(*(user->mi), {i,j});
        f[j][i] += user->trPtr->diffusion::D * 
                   user->trPtr->diffusion::Flux(*(user->mi), {i,j});
    }}

    //! Restore array to local vectors.
    DMDAVecRestoreArray(dmu, F, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}

/**
 * Explicit advection.
 * Using time stepping object provided by Petsc.
 */
PetscErrorCode ExplicitAdvection(TS ts, PetscReal time, Vec U, Vec F, void* ctx){

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

    // Update smoothness indicators for all reconstruction levels altogether
    user->trPtr->UpdateSmoothnessIndic(*(user->mi));

    // Update non linear weights based on updated smoothness indicators
    user->trPtr->advection::UpdateNonLinearWgts(*(user->mi));

    // Update flux defined on each edge
    user->trPtr->advection::UpdateEdgeFlux(*(user->mi));

    //! Loop through computational domain
    for (int j=user->mi->MPIlocalCellStart[1]; j<user->mi->MPIlocalCellStart[1] + user->mi->MPIlocalCellSize[1]; j++){
    for (int i=user->mi->MPIlocalCellStart[0]; i<user->mi->MPIlocalCellStart[0] + user->mi->MPIlocalCellSize[0]; i++){
        f[j][i] = -1.0*user->trPtr->advection::Flux(*(user->mi), {i,j});
    }}

    //! Restore array to local vectors.
    DMDAVecRestoreArray(dmu, F, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}

/**
 * Explicit diffusion.
 * Using time stepping object provided by Petsc.
 */
PetscErrorCode ExplicitDiffusion(TS ts, PetscReal time, Vec U, Vec F, void* ctx){

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

    // Update smoothness indicators for all reconstruction levels altogether
    user->trPtr->UpdateSmoothnessIndic(*(user->mi));

    // Update non linear weights based on updated smoothness indicators
    user->trPtr->diffusion::UpdateNonLinearWgts(*(user->mi));

    // Update flux defined on each edge
    user->trPtr->diffusion::UpdateEdgeFlux(*(user->mi));

    //! Loop through computational domain
    for (int j=user->mi->MPIlocalCellStart[1]; j<user->mi->MPIlocalCellStart[1] + user->mi->MPIlocalCellSize[1]; j++){
    for (int i=user->mi->MPIlocalCellStart[0]; i<user->mi->MPIlocalCellStart[0] + user->mi->MPIlocalCellSize[0]; i++){
        f[j][i] = user->trPtr->diffusion::D * 
                  user->trPtr->diffusion::Flux(*(user->mi), {i,j});
    }}

    //! Restore array to local vectors.
    DMDAVecRestoreArray(dmu, F, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}

/**
 * Compute jacobian requiared for implicit time stepping
 * Advection-diffusion equation
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

    // Update edge flux first
    user->trPtr->advection::UpdateEdgeFluxDerivative(*(user->mi));
    user->trPtr->diffusion::UpdateEdgeFluxDerivative(*(user->mi));

    //! Get MPI local part of Jacobian matrix
    int rstart, rend;
    MatGetOwnershipRange(J, &rstart, &rend);

    unordered_map<int,double> deriv;

    for (int row = rstart; row<rend; row++){
        indice global = Bend(*(user->mi), row);

        double opp = -1.0;

        // add advection flux to deriv
        unordered_map_arithmetic(deriv, 
                                 user->trPtr->advection::derivFlux(*(user->mi),global),
                                 std::minus<double>());

        // add diffusion flux to deriv
        unordered_map_arithmetic(deriv,
                                 user->trPtr->diffusion::derivFlux(*(user->mi),global),
                                 std::plus<double>(),
                                 user->trPtr->diffusion::D,
                                 std::multiplies<double>());

        for (auto & dVal: deriv){
            ierr = MatSetValue(J,row,dVal.first,dVal.second,INSERT_VALUES);CHKERRQ(ierr);
        }

        deriv.clear();

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

/**
 * Compute jacobian requiared for implicit time stepping
 * Pure advection equation
 */
PetscErrorCode FormJacobianAdvection(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void* ctx){
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

    // Update edge flux first
    user->trPtr->advection::UpdateEdgeFluxDerivative(*(user->mi));

    //! Get MPI local part of Jacobian matrix
    int rstart, rend;
    MatGetOwnershipRange(J, &rstart, &rend);

    unordered_map<int,double> deriv;

    for (int row = rstart; row<rend; row++){
        indice global = Bend(*(user->mi), row);

        // advection flux
        deriv = user->trPtr->advection::derivFlux(*(user->mi), global);

        for (auto & dVal: deriv){
            ierr = MatSetValue(J,row,dVal.first,-1*dVal.second,INSERT_VALUES);CHKERRQ(ierr);
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

/**
 * Compute jacobian requiared for implicit time stepping
 * Pure diffusion equation
 */
PetscErrorCode FormJacobianDiffusion(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void* ctx){
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

    // Update edge flux first
    user->trPtr->diffusion::UpdateEdgeFluxDerivative(*(user->mi));

    //! Get MPI local part of Jacobian matrix
    int rstart, rend;
    MatGetOwnershipRange(J, &rstart, &rend);

    unordered_map<int,double> deriv;

    for (int row = rstart; row<rend; row++){
        indice global = Bend(*(user->mi), row);

        // Diffusion flux
        deriv = user->trPtr->diffusion::derivFlux(*(user->mi), global);

        for (auto & dVal: deriv){
            ierr = MatSetValue(J,row,dVal.first,user->trPtr->diffusion::D * dVal.second,INSERT_VALUES);CHKERRQ(ierr);
        }

    }

    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (J != Jp){
        ierr = MatAssemblyBegin(Jp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Jp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

    // Output Jacobian
    std::string tmp = std::to_string(time);

    DrawMat(J,tmp.c_str());

    //! Restore array to local vectors
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}

/**
 * Explicit advection-diffusion.
 * Using time stepping object provided by Petsc.
 */
PetscErrorCode ExplicitFull(TS ts, PetscReal time, Vec U, Vec F, void* ctx){

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

    // Update smoothness indicators for all reconstruction levels altogether
    user->trPtr->UpdateSmoothnessIndicAndDerivative(*(user->mi));

    // Update non linear weights based on updated smoothness indicators
    user->trPtr->advection::UpdateNonLinearWgtsAndDerivative(*(user->mi));
    user->trPtr->diffusion::UpdateNonLinearWgtsAndDerivative(*(user->mi));

    // Update flux defined on each edge
    user->trPtr->advection::UpdateEdgeFlux(*(user->mi));
    user->trPtr->diffusion::UpdateEdgeFlux(*(user->mi));

    //! Loop through computational domain
    for (int j=user->mi->MPIlocalCellStart[1]; j<user->mi->MPIlocalCellStart[1] + user->mi->MPIlocalCellSize[1]; j++){
    for (int i=user->mi->MPIlocalCellStart[0]; i<user->mi->MPIlocalCellStart[0] + user->mi->MPIlocalCellSize[0]; i++){
        f[j][i] = -1.0*user->trPtr->advection::Flux(*(user->mi), {i,j});
        f[j][i] += user->trPtr->diffusion::D * 
                   user->trPtr->diffusion::Flux(*(user->mi), {i,j});
    }}

    //! Restore array to local vectors.
    DMDAVecRestoreArray(dmu, F, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}

/**
 * Compute jacobian requiared for implicit time stepping
 * Advection-diffusion equation.
 * Full differentiation including non linear weights.
 */
PetscErrorCode FormJacobianFull(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void* ctx){
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

    // Update edge flux first
    user->trPtr->advection::UpdateEdgeFluxDerivativeFull(*(user->mi));
    user->trPtr->diffusion::UpdateEdgeFluxDerivativeFull(*(user->mi));

    //! Get MPI local part of Jacobian matrix
    int rstart, rend;
    MatGetOwnershipRange(J, &rstart, &rend);

    derivative deriv;

    for (int row = rstart; row<rend; row++){
        indice global = Bend(*(user->mi), row);

        // advection flux
        deriv = user->trPtr->advection::derivFlux(*(user->mi), global);

        for (auto & dVal: deriv){
            ierr = MatSetValue(J,row,dVal.first,-1*dVal.second,ADD_VALUES);CHKERRQ(ierr);
        }

        // Diffusion flux
        deriv = user->trPtr->diffusion::derivFlux(*(user->mi), global);

        for (auto & dVal: deriv){
            ierr = MatSetValue(J,row,dVal.first,user->trPtr->diffusion::D * dVal.second,ADD_VALUES);CHKERRQ(ierr);
        }

    }

    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (J != Jp){
        ierr = MatAssemblyBegin(Jp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Jp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

    // Output Jacobian
    std::string tmp = std::to_string(time);

    DrawMat(J,tmp.c_str());

    //! Restore array to local vectors
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    PetscFunctionReturn(0);
}
