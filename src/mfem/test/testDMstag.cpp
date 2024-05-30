#include <iostream>
#include <fstream>
#include <petsc.h>

extern "C"{
#include "mesh.h"
}

int main(int argc, char ** argv){

    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n"));

    int M = 2, N = 2;
    PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);
    PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);

    DM test;

    DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, 
                   M, N, PETSC_DECIDE, PETSC_DECIDE, 2, 1, 0, DMSTAG_STENCIL_BOX,0,
                   NULL, NULL,&test);
    DMSetUp(test); 

    DMView(test, PETSC_VIEWER_STDOUT_WORLD);

    // ======================
    Vec   fullmesh;
    DM    dm;

    const int stencilWidth = 2;

    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidth, NULL, NULL, &dm);
    DMSetFromOptions(dm); 
    DMSetUp(dm); 
    DMCreateGlobalVector(dm, &fullmesh); 

    MeshParam mp;
    mp.xstart = -1;
    mp.ystart = -1;
    mp.L = 2;
    mp.H = 2;

    CreateFullMesh(dm, &fullmesh, &mp);

    VecDestroy(&fullmesh);
    DMDestroy(&dm);
    DMDestroy(&test);

    PetscFinalize();

    return 0;
}
