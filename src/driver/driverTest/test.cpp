#include "driver.h"

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

    double L = 0.2, H = 1;
    double xstart = -0.1, ystart = -1.0001;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    int stencilWidthMesh = 5;
    int stencilWidthU = 3;

    int physicsScale = 1;
    PetscCall(PetscOptionsGetInt(NULL,NULL, "-scale", &physicsScale, NULL));

    int meshType = 0; 
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshType,NULL));

    int maxIter = 1; 
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-maxIter", &maxIter, NULL));        
    double tolUzawa = 10e-7; 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tol", &tolUzawa, NULL)); 

    double Tmax = 0.009; // Stop at the first step 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tmax", &Tmax, NULL)); 

    // ==============================================================================

    Driver * driver = new Driver();

    driver->CreatePhase();

    driver->CreateMesh(M, N, L, H, xstart, ystart, 
                       stencilWidthMesh, stencilWidthU,
                       physicsScale, meshType);

    driver->InitTransport(InitHD, InitCD);

    //cout << InitHD({0,-1},{driver->myPhase->pp->l0,0.0}) << endl;

    driver->PrepareDefaultTransport();

    driver->PrepareFlow();

    driver->SolveFlow(maxIter, tolUzawa);

    driver->CreateScatterVec();

    // Print initial condition
    driver->PrintGrid();
    driver->PrintLithoPressure();

    driver->Tmax = Tmax;
    driver->maxIter = maxIter;
    driver->dt = 0.01;

    driver->RK();

/*
    // Time stepping
    ctx_driver ctx;
    ctx.driver = driver;
    ctx.dt = 0.01;
    ctx.maxIter  = maxIter;
    ctx.tolUzawa = tolUzawa;

    TS ts;
    PetscCall(TSCreate(PETSC_COMM_WORLD, &ts)); 
    TSSetProblemType(ts, TS_NONLINEAR);
    //TSSetMaxTime(ts, 0.2);
    TSSetMaxTime(ts, 1.0);

    TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
    TSSetDM(ts, driver->dmu);
    TSSetTimeStep(ts, ctx.dt);

    Vec U;
    Vec array[2];
    array[0] = driver->globalCD;
    array[1] = driver->globalHD;
    PetscCall(VecCreateNest(PETSC_COMM_WORLD,2,NULL,array,&U));

    TSSetSolution(ts, U);

    TSSetRHSFunction(ts, NULL, Explicit, &ctx);

    TSSetType(ts, TSEULER);

    TSSetUp(ts);

    TSSolve(ts, U);
*/

    // =============== Print functions ==============================================

    driver->eventCount++;

    driver->PrintPhaseEvent();
    driver->PrintFlowEvent();
    driver->PrintPressureEvent();
    driver->PrintHDEvent();
    driver->PrintCDEvent();

    //driver->clean();

    PetscFinalize();

    return 0;
}
