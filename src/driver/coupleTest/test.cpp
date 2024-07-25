#include <petsc.h>
#include <iostream>
#include "integral.h"
#include "phase.h"
#include "param.h"

// MFEM parameter header file
//#include "myFunc.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

int main(int argc, char **argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Initialize mesh ==========================================================
    int M=3, N=3;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL));

    PhysProperty * pp = (PhysProperty *)malloc(sizeof(PhysProperty));

    AssignPhyProperties(pp);

    // Physical domain
    double physscale = pp->L0/pp->l0;
    double L = 2*physscale, H = 1*physscale;
    double xstart = -1*physscale, ystart = -1.0001*physscale;

    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    DM dmMesh;

    int stencilWidthMesh = 5;

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidthMesh, NULL, NULL, &dmMesh));
    PetscCall(DMSetFromOptions(dmMesh));              
    PetscCall(DMSetUp(dmMesh));

    Vec globalmesh;

    PetscCall(DMCreateGlobalVector(dmMesh, &globalmesh));

    // If running test on the single stencil for convergence study
    int singleStencilTest = 0;
    double scale = 1;
    PetscCall(PetscOptionsGetInt(NULL,NULL, "-single", &singleStencilTest, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL, "-scale", &scale, NULL));

    if (singleStencilTest){
        L = L/scale;
        H = H/scale;
        xstart = -L/2.0;
        ystart = -H/2.0;
    }

    // Define Intermediate data structure
    MeshParam mp;
    mp.xstart = xstart;
    mp.ystart = ystart;
    mp.L = L;
    mp.H = H;

    int meshtype = 0; 

    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshtype,NULL));
    switch(meshtype){
        case 0: CreateFullMesh(dmMesh, &globalmesh, &mp); break;
        case 1: LogicRectMesh(dmMesh, &globalmesh, &mp);  break;
        case 2: RefineMesh(dmMesh, &globalmesh, &mp);
        //case 2: TestControlMeshSecond(dmCell,L,H); break;
        //case 3: TestControlMeshThird(dmCell,L,H);  break;
    }

    // Print mesh or not
    int printmesh=0;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-printmesh",&printmesh,NULL));
    if(printmesh){ 
        VecView(globalmesh, PETSC_VIEWER_STDOUT_WORLD);
        PrintFullMesh(dmMesh, &globalmesh);
    }

    // Setting initial solution
    DM dmu; 
    int stencilWidthU = 3;

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 1, 
    stencilWidthU, NULL, NULL, &dmu));
    PetscCall(DMSetFromOptions(dmu));              
    PetscCall(DMSetUp(dmu));     

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n"));
    // ==========================================================================

    // Create initial (C,H) distribution pair
    // Calculate volume fraction of fluid (porosity)
    PorosityOut(xstart, ystart, L, H, 20, pp);      


    // Solve for velocity

    // Transport


    // Finialize the program ====================================================
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    PetscFinalize();

    return 0;
}
