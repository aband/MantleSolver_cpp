#include "driver.h"

int Driver::CreatePhase(){

    myPhase = new Phase();

    myPhase->pp = (PhysProperty *)malloc(sizeof(PhysProperty));

    AssignPhyProperties(myPhase->pp);

    myPhase->pPtr = new EUTECTIC::phase();

    return 0;
}

int Driver::CreateDMs(const int& M, const int& N,
                      double L, double H, 
							 double xstart, double ystart,
                      const int& stencilWidthMesh, 
							 const int& stencilWidthU,
					       const bool& physicsScale){

    if (physicalScale){
        double physscale = myPhase->pp->L0/myPhase->pp->l0;
        L = L*physscale;
        H = H*physscale;
        xstart = xstart*physscale, 
        ystart = ystart*physscale;
	 }

    // Create dmMesh
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidthMesh, NULL, NULL, &dmMesh));
    PetscCall(DMSetFromOptions(dmMesh));              
    PetscCall(DMSetUp(dmMesh));

    // Create dmU
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 1, 
    stencilWidthU, NULL, NULL, &dmu));
    PetscCall(DMSetFromOptions(dmu));              
    PetscCall(DMSetUp(dmu));     

    // Create MeshParam object (historical object one time use only)
    mp_.xstart = xstart;
    mp_.ystart = ystart;
    mp_.L = L;
    mp_.H = H;

    return 0;
}
