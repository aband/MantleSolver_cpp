#include <iostream>
#include <petsc.h>
#include "integral.h"

// From MLWENO
#include "stencil.h"
#include "util.h"
#include "input.h"
#include "reconstruction.h"

// From tranport
#include "transport.h"
#include "timestepping.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

#include <chrono>

using namespace std;

double InitialValue(vertex& point, const vector<double>& param){
	 //if (point[0]<-1.0/param[0]){
//		  return point[0]*point[0]+point[1]*point[1];
//	     return sin(point[0]*3.0)+cos(point[1]/2.0) + point[0]*(point[1]+1);
//	 } else {
//		  return point[0]*point[0]*point[1]*point[1] + 1.0;
//	     return sin(point[0]*3.0)+cos(point[1]/2.0) + point[0]*(point[1]+1) + 1;
//	 }

    //return sin(point[0]*3.0+0.5)+cos(point[1]/2.0-0.2) + pow(point[0]+0.1,3)*(point[1]+1);
    //return point[0]*point[0] + point[1]*point[1];
    //return point[0] + point[1];

    // Initial value for sine wave 2D Burger's equation
    return pow(sin(M_PI*(point[0]+1)/2),2)*pow(sin(M_PI*(point[1]+1)/2),2);

}

PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal t, Vec U, void *ctx){

    PetscFunctionBeginUser;

    SNES snes;
    TSGetSNES(ts, &snes);

    int it;
    SNESGetIterationNumber(snes, &it);

    cout << "T = " << t << ". Newton iteration number: " << it << endl;

    PetscFunctionReturn(0);
}

int main(int argc, char **argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size);CHKERRQ(ierr);

    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

// ==========================================================================================================================

    // Start testing mesh function
    // Initializing problem size with 3X3
    int M = 5, N = 5;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    // Create data management object
    DM    dm;
    Vec   fullmesh;
    const int stencilWidth = 5;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidth, NULL, NULL, &dm);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dm);               CHKERRQ(ierr);
    ierr = DMSetUp(dm);                        CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm, &fullmesh);CHKERRQ(ierr); 

    double L = 2.0, H = 2.0;
    double xstart = -1.0, ystart = -1.0;
    ierr = PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL); CHKERRQ(ierr);

    MeshParam mp;
    mp.xstart = xstart;
    mp.ystart = ystart;
    mp.L = L;
    mp.H = H;

    // Uniform or distorted mesh
    int meshtype=0;
    ierr = PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshtype,NULL);CHKERRQ(ierr);
    switch(meshtype){
        case 0: CreateFullMesh(dm, &fullmesh, &mp); break;
        case 1: LogicRectMesh(dm, &fullmesh, &mp);  break;
        case 2: RefineMesh(dm, &fullmesh, &mp); break;
        //case 2: TestControlMeshSecond(dmCell,L,H); break;
        //case 3: TestControlMeshThird(dmCell,L,H);  break;
    }

    int printmesh=0;
    ierr = PetscOptionsGetInt(NULL,NULL,"-printmesh",&printmesh,NULL);CHKERRQ(ierr);
    if(printmesh){ 
        VecView(fullmesh, PETSC_VIEWER_STDOUT_WORLD);
        PrintFullMesh(dm, &fullmesh);
    }

    cout << "Mesh Created. To check full mesh, rerun with -printmesh 1 " << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

// ==========================================================================================================================

    // Contain defined mesh in vector container.
    // and verify it.
    vector< valarray<double> > mesh;
    
    ReadMeshPortion(dm, &fullmesh, mesh);

    //cout << "Converted c array of local mesh into vector container c++ " << endl;
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

// ====================================================================================================================================

    DM dmu;

    int cell_ghost = 3;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 1, cell_ghost, NULL, NULL, &dmu);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dmu); CHKERRQ(ierr);
    ierr = DMSetUp(dmu);          CHKERRQ(ierr);

    Vec globalu;
    ierr = DMCreateGlobalVector(dmu,&globalu);CHKERRQ(ierr);

    // Initialize with oblique data for Burgers equation 
    //ObliqueBurgers(dm,dmu,&fullmesh,&globalu,Initial_Condition);
    SimpleInitialValue(dm,dmu,&fullmesh,&globalu,InitialValue);

    Vec localu; 
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, localu);

    // It can be changed later to not be double
    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

// ====================================================================================================================================

    // Create MeshInfo object
    MeshInfo mi; 

    mi.lmesh = mesh;
    mi.localVals = lu;

    // Assign meshinfo after mesh added to meshinfo
    AssignValuesMeshInfo(mi,dm,dmu); 

// ========================================================================================================================================

    // Ouptut of initial value
    char * filename = (char *)"initial.txt";

    PlainOutput(dmu, &globalu, filename);
    PlainMeshOutput(dm, &fullmesh);

// ========================================================================================================================================

    //! Create a transport object
    transport* trPtr = new transport();

    /**
     * Initialize transport object
     * No need to give a precious description of reconstruction methods.
     * Detailed reconstruction method will be added separately later.
     */
    trPtr->AddLevel(mi,1,1);
    trPtr->AddLevel(mi,2,2);
    trPtr->AddLevel(mi,3,3);
    trPtr->AddLevel(mi,4,5);
    trPtr->AddLevel(mi,5,4);

    trPtr->AssignReconstMethod(trPtr->advection::reconstMethods,"(1,1)",{{0,0}});
    trPtr->AssignReconstMethod(trPtr->advection::reconstMethods,"(2,2)",{{-1,0},{0,0},{-1,-1},{0,-1}});
    trPtr->AssignReconstMethod(trPtr->advection::reconstMethods,"(3,3)",{{-1,-1}});

    trPtr->AssignReconstMethod(trPtr->diffusion::reconstMethodsVert,"(1,1)",{{0,0}});
    trPtr->AssignReconstMethod(trPtr->diffusion::reconstMethodsVert,"(2,2)",{{-1,0},{0,0},{-1,-1},{0,-1}});
    //trPtr->AssignReconstMethod(trPtr->diffusion::reconstMethodsVert,"(3,3)",{{-1,0},{-2,0},{-2,-2},{-1,-2}});
    //trPtr->AssignReconstMethod(trPtr->diffusion::reconstMethodsVert,"(4,5)",{{-2,-2}});

    trPtr->AssignReconstMethod(trPtr->diffusion::reconstMethodsHori,"(1,1)",{{0,0}});
    trPtr->AssignReconstMethod(trPtr->diffusion::reconstMethodsHori,"(2,2)",{{-1,0},{0,0},{-1,-1},{0,-1}});
    //trPtr->AssignReconstMethod(trPtr->diffusion::reconstMethodsHori,"(3,3)",{{0,-1},{-2,-1},{-2,-2},{0,-2}});
    //trPtr->AssignReconstMethod(trPtr->diffusion::reconstMethodsHori,"(5,4)",{{-2,-2}});

    trPtr->CreateMLWENO(mi);

// Test and print ========================================================
    //trPtr->UpdateSmoothnessIndic(mi);

    //trPtr->advection::UpdateNonLinearWgts(mi);

    //trPtr->advection::GetInfo(mi);

    //trPtr->diffusion::UpdateNonLinearWgts(mi);

    //trPtr->diffusion::GetInfo(mi);
// =======================================================================

    //! Create ctx for time stepping
    Ctx ctx;
    ctx.trPtr = trPtr;
    ctx.mi    = &mi;
    ctx.dmu   = dmu;

    /**
     * Time stepping.
     */

    TS ts;
    SNES snes;
    KSP ksp;
    PC  pc;

    double Tmax = 0.01;
    double dt = 0.01;

    //! Get time variables input from terminal line
    ierr = PetscOptionsGetReal(NULL,NULL,"-Tmax",&Tmax,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);    CHKERRQ(ierr);
 
    //SNESCreate(PETSC_COMM_WORLD, &snes);
    //KSPCreate(PETSC_COMM_WORLD, &ksp);

    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetProblemType(ts, TS_NONLINEAR);

    TSSetMaxTime(ts, Tmax);
    TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
    TSSetDM(ts,dmu);

    TSSetTimeStep(ts, dt);
    TSSetSolution(ts,globalu);

    SNESLineSearch linesearch;

    //! Change preconditioner
    TSGetSNES(ts, &snes);
    SNESGetKSP(snes, &ksp);
    SNESGetLineSearch(snes, &linesearch);
    SNESLineSearchSetTolerances(linesearch, 0, 1e8, 1e-8, 1e-14, 1e-8, 30);
    KSPGetPC(ksp, &pc);
    KSPSetTolerances(ksp, 1e-8,1e-13,1000,30);
    PCSetType(pc, PCJACOBI);
    PCSetFromOptions(pc);

    //TSSetRHSFunction(ts, globalu, Explicit, &ctx);
    TSSetRHSFunction(ts, NULL, Explicit, &ctx);

    // ===================================================================
    //! Explicit
    //! Forward Euler
    TSSetType(ts, TSEULER);

    //! SSP
    //TSSetType(ts, TSSSP);
    //TSSSPSetType(ts, TSSSPRKS2);

    //TSRKSetType(ts, TSRK3);

    // ===================================================================
    //! Implicit
/*
    Mat J;
    MatCreate(PETSC_COMM_WORLD, &J);
    MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, M*N, M*N);
    MatSetUp(J);
    TSSetType(ts, TSBEULER);
    TSSetRHSJacobian(ts, J, J, FormJacobian, &ctx);
*/

    //TSSetTolerances(ts,1e-3,NULL,1e-3,NULL);
    //TSSetMaxSNESFailures(ts, 50);

    TSMonitorSet(ts, Monitor, &ctx, NULL);

    TSSetFromOptions(ts);
    TSSetUp(ts);

    cout << "Time stepping begins .. .. .. " << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    auto start = std::chrono::system_clock::now();
    TSSolve(ts,globalu);
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;

    cout << "Elapsed time: " << elapsed_seconds.count() << endl;

    TSView(ts,PETSC_VIEWER_STDOUT_WORLD);

// ====================================================================================================================================
    // Ouptut of final result
    filename = (char *)"final.txt";

    PlainOutput(dmu, &globalu, filename);

    delete trPtr;

// ====================================================================================================================================
    // Clear used objects
    DMDAVecRestoreArray(dmu,localu,&lu);
    DMRestoreLocalVector(dmu, &localu); 

    VecDestroy(&fullmesh);
    VecDestroy(&globalu);
    DMDestroy(&dm);
    DMDestroy(&dmu);

    return 0;
}
