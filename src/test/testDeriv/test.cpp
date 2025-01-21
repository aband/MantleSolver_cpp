// Test calculation of derivatives of WENO reconstruction

#include "driver.h"

int main(int argc, char **argv){

    PetscMPIInt size, rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    int M = 3, N = 3; 

    PetscCall(PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL));

    double L = 1, H = 1;
    double xstart = 0.0, ystart = 0.0;
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

    // =====================================================================
    Driver * driver = new Driver();

    driver->getDomainSize(L, H);

    driver->CreateMesh(M, N, L, H, xstart, ystart, 
                       stencilWidthMesh, stencilWidthU,
                       physicsScale, meshType);

    driver->InitTransport(InitCD, "test", false);

    //driver->PrintMesh("test");

    driver->AddLevels(1);
    driver->AddLevels(2);
    driver->AddLevels(3);

    driver->PrepareTransport("test");



    return 1;
}
