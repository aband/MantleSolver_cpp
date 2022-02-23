#include <petsc.h>
#include "../include/mesh.h"

PetscErrorCode CreateFullMesh(DM dm, Vec *fullmesh, MeshParam * mp){

    PetscErrorCode ierr;
    Vec            fmesh, lmesh;
    PetscInt       xs,ys,xm,ym,M,N;
    Point          **localmesh;
    PetscInt       stencilwidth;
    PetscFunctionBeginUser;

    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                                CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Unpack parameters
    double L = mp->L;
    double H = mp->H;
    double xstart = mp->xstart;
    double ystart = mp->ystart;

    double hx = L/(double)M;
    double hy = H/(double)N;

    // Define rectangular mesh
    ierr = DMGetLocalVector(dm, &lmesh);                         CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dm, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr); 
    ierr = DMGlobalToLocalEnd(dm, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dm, lmesh, &localmesh);               CHKERRQ(ierr);

    for (int j=ys-stencilwidth; j<ys+ym+stencilwidth; j++){
    for (int i=xs-stencilwidth; i<xs+xm+stencilwidth; i++){
        localmesh[j][i].p[0] = xstart+i*hx;
        localmesh[j][i].p[1] = ystart+j*hy;
    }}

    ierr = DMDAVecRestoreArray(dm, lmesh, &localmesh); CHKERRQ(ierr);

    // Update global vector with ghost region
    ierr = DMLocalToGlobalBegin(dm, lmesh, ADD_VALUES, fmesh);CHKERRQ(ierr); 
    ierr = DMLocalToGlobalEnd(dm, lmesh, ADD_VALUES, fmesh);  CHKERRQ(ierr);

    // Global to Local again
    ierr = DMGlobalToLocalBegin(dm, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr); 
    ierr = DMGlobalToLocalEnd(dm, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dm, lmesh, &localmesh);               CHKERRQ(ierr);

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+xm; i++){
        localmesh[j][i].p[0] = xstart+i*hx;
        localmesh[j][i].p[1] = ystart+j*hy;
    }}

    ierr = DMDAVecRestoreArray(dm, lmesh, &localmesh); CHKERRQ(ierr);

    // Update global vector without ghost region
    ierr = DMLocalToGlobalBegin(dm, lmesh, INSERT_VALUES, fmesh);CHKERRQ(ierr); 
    ierr = DMLocalToGlobalEnd(dm, lmesh, INSERT_VALUES, fmesh);  CHKERRQ(ierr);

    ierr = DMRestoreLocalVector(dm, &lmesh);           CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode LogicRectMesh(DM dm, Vec *fullmesh, MeshParam * mp){

    PetscErrorCode ierr;
    Vec            fmesh, lmesh;
    PetscInt       xs,ys,xm,ym,M,N;
    Point          **localmesh;
    PetscInt       stencilwidth;
    PetscFunctionBeginUser;

    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                                CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Unpack parameters
    double L = mp->L;
    double H = mp->H;
    double xstart = mp->xstart;
    double ystart = mp->ystart;

    double hx = L/(double)M;
    double hy = H/(double)N;

    // Define rectangular mesh
    ierr = DMGetLocalVector(dm, &lmesh);                         CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dm, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr); 
    ierr = DMGlobalToLocalEnd(dm, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dm, lmesh, &localmesh);               CHKERRQ(ierr);

    for (int j=ys-stencilwidth; j<ys+ym+stencilwidth; j++){
    for (int i=xs-stencilwidth; i<xs+xm+stencilwidth; i++){
        localmesh[j][i].p[0] = xstart+i*hx;
        localmesh[j][i].p[1] = ystart+j*hy;
    }}

    ierr = DMDAVecRestoreArray(dm, lmesh, &localmesh); CHKERRQ(ierr);

    // Update global vector with ghost region
    ierr = DMLocalToGlobalBegin(dm, lmesh, ADD_VALUES, fmesh);CHKERRQ(ierr); 
    ierr = DMLocalToGlobalEnd(dm, lmesh, ADD_VALUES, fmesh);  CHKERRQ(ierr);

    // Global to Local again
    ierr = DMGlobalToLocalBegin(dm, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr); 
    ierr = DMGlobalToLocalEnd(dm, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dm, lmesh, &localmesh);               CHKERRQ(ierr);

    PetscRandom    rndx,rndy;
    PetscScalar    value;

    // Get current time in seconds as random seed
    time_t seconds;
    time(&seconds);

    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rndx);CHKERRQ(ierr);
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rndy);CHKERRQ(ierr);

    ierr = PetscRandomSetInterval(rndx,-0.3*hx,0.3*hx)    ;CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(rndx,(unsigned long)seconds);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rndx)                ;CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(rndy,-0.3*hy,0.3*hy)    ;CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(rndy,(unsigned long)seconds);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rndy)                ;CHKERRQ(ierr);

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+xm; i++){
        if (i*j != 0){
            ierr = PetscRandomGetValue(rndy, &value);CHKERRQ(ierr);
            localmesh[j][i].p[0] = xstart+i*hx + (double)(PetscRealPart(value));
            ierr = PetscRandomGetValue(rndx, &value);CHKERRQ(ierr);
            localmesh[j][i].p[1] = ystart+j*hy + (double)(PetscRealPart(value));
        }
    }}

    ierr = DMDAVecRestoreArray(dm, lmesh, &localmesh); CHKERRQ(ierr);

    // Update global vector without ghost region
    ierr = DMLocalToGlobalBegin(dm, lmesh, INSERT_VALUES, fmesh);CHKERRQ(ierr); 
    ierr = DMLocalToGlobalEnd(dm, lmesh, INSERT_VALUES, fmesh);  CHKERRQ(ierr);

    ierr = DMRestoreLocalVector(dm, &lmesh);           CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
