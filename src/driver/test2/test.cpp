#include "driver.h"
#include "print.h"

int main(int argc, char **argv){

    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Input mesh parameter =========================================================
    int M=4, N=20;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL));

    double L = 0.1, H = 0.4;
    double xstart = -0.5*L, ystart = -1.0001*H;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    int stencilWidthMesh = 5;
    int stencilWidthU = 3;

    int physicsScale = 0;
    PetscCall(PetscOptionsGetInt(NULL,NULL, "-scale", &physicsScale, NULL));

    int meshType = 0; 
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshType,NULL));

    int maxIter = 15; 
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-maxIter", &maxIter, NULL));        
    double tolUzawa = 10e-15; 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tol", &tolUzawa, NULL)); 

    double Tmax = 20; // Stop at the first step 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tmax", &Tmax, NULL)); 

    double dt = 1;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL));

    int showPhase = 0;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-showphase", &showPhase, NULL)); 

    int withUnit = 0;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-unit", &withUnit, NULL));

    // ==============================================================================

    Driver * driver = new Driver();

    driver->withUnit = withUnit;
    driver->CreatePhase();
    driver->ShowPhase();
    driver->CreateMesh(M, N, L, H, xstart, ystart, 
                       stencilWidthMesh, stencilWidthU,
                       physicsScale, meshType);

    /**!
     * Initialize global cell averaged value vectors.
     * Initialize multi level reconstruction objects
     */
    driver->PrepareTransport(InitHD, InitCD);

    printCellCenterGrid(driver->mi);
    printCellAve(1, &driver->globalHD, driver->mi, "HD");
    printCellAve(1, &driver->globalCD, driver->mi, "CD");

    /**!
     * Create boundary reference arrays
     * Allocate memory space for solutions vectors
     */
    driver->PrepareFlow();

    /**!
     * One step computation for flow problem
     */
    double h0 = sqrt((L*H)/(double)(M*N));

    // Solve for initial velocity
    Vec localHD, localCD;
    double ** lHD;
    double ** lCD;

    PetscCall(DMGetLocalVector(driver->dmu, &localHD)); 

    PetscCall(DMGlobalToLocalBegin(driver->dmu, driver->globalHD, INSERT_VALUES, localHD));
    PetscCall(DMGlobalToLocalEnd(driver->dmu, driver->globalHD, INSERT_VALUES, localHD));

    PetscCall(DMDAVecGetArray(driver->dmu, localHD, &lHD));

    PetscCall(DMGetLocalVector(driver->dmu, &localCD)); 

    PetscCall(DMGlobalToLocalBegin(driver->dmu, driver->globalCD, INSERT_VALUES, localCD));
    PetscCall(DMGlobalToLocalEnd(driver->dmu, driver->globalCD, INSERT_VALUES, localCD));

    PetscCall(DMDAVecGetArray(driver->dmu, localCD, &lCD));

    driver->ml.updatesigma(lHD);
    Tensor<weights> allwgtsHD;
    driver->advection.computeWgts(driver->ml, driver->mi, h0, allwgtsHD, location);

    driver->ml.updatesigma(lCD);
    Tensor<weights> allwgtsCD;
    driver->advection.computeWgts(driver->ml, driver->mi, h0, allwgtsCD, location);

    driver->V0 = driver->myPhase->pp->V0 / driver->myPhase->pPtr->u0;
    cout << driver->V0 << endl;

    driver->SolveFlow(maxIter, tolUzawa, allwgtsHD, lHD, allwgtsCD, lCD);

    driver->PrintFlowEvent(1);
    //driver->PrintFlowEventTransform(1);

    //driver->PrintPhaseEvent(1);

    driver->PrintEffVel(1, 2, allwgtsHD, lHD, allwgtsCD, lCD);

    DMDAVecRestoreArray(driver->dmu,localHD,&lHD);
    DMRestoreLocalVector(driver->dmu, &localHD); 
    DMDAVecRestoreArray(driver->dmu,localCD,&lCD);
    DMRestoreLocalVector(driver->dmu, &localCD); 

    driver->PrintPressureSerialApprox(1);

    /**!
     * Actual time stepping.
     */
//    driver->RK(dt, Tmax, maxIter, tolUzawa);

    VecDestroy(&driver->globalmesh);
    VecDestroy(&driver->globalHD);
    VecDestroy(&driver->globalCD);
    DMDestroy(&driver->dmu);
    DMDestroy(&driver->dmMesh);

    PetscFinalize();

    return 0;
}
