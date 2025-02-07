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

    int M = 3, N = 3;
    PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);
    PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);

    // ======================

    int stencilwidth = 0;

    DM    dmbr;
    DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, M, N, PETSC_DECIDE, PETSC_DECIDE, 2, 1, 0, 
                   DMSTAG_STENCIL_BOX, stencilwidth, NULL, NULL, &dmbr);
    DMSetFromOptions(dmbr); 
    DMSetUp(dmbr); 

    DMView(dmbr, PETSC_VIEWER_STDOUT_WORLD);

    Vec testbr;
    DMCreateGlobalVector(dmbr, &testbr);

    Mat Abr;
    DMCreateMatrix(dmbr, &Abr);

    // ==================================================

    DM dmbdm;
    DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, M, N, PETSC_DECIDE, PETSC_DECIDE, 0, 2, 0, 
                   DMSTAG_STENCIL_BOX, stencilwidth, NULL, NULL, &dmbdm);
    DMSetFromOptions(dmbdm); 
    DMSetUp(dmbdm); 

    DMView(dmbdm, PETSC_VIEWER_STDOUT_WORLD);

    Vec testbdm;
    DMCreateGlobalVector(dmbdm, &testbdm);

    Mat Abdm;
    DMCreateMatrix(dmbdm, &Abdm);

    DMDestroy(&dmbr);
    DMDestroy(&dmbdm);

    PetscFinalize();

    return 0;
}
