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

    PetscInt startx, starty, nx, ny;

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

    PetscCall(DMStagGetCorners(dmbr, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
    PetscCall(DMStagGetGlobalSizes(dmbr, &M, &N, NULL));

    DMStagStencil sten;
    sten.i = 1;
    sten.j = 1;
    sten.loc = DMSTAG_UP;
    sten.c = 0;

    double val = 1.0;

    DMStagMatSetValuesStencil(dmbr, Abr, 1, &sten, 1, &sten, &val, INSERT_VALUES);

/*
    for (int j=starty; j<starty+ny; j++){
    for (int i=startx; i<startx+nx; i++){

        int dof[12];


    }}
*/

    MatAssemblyBegin(Abr, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Abr, MAT_FINAL_ASSEMBLY);

    MatView(Abr, PETSC_VIEWER_STDOUT_WORLD);

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
