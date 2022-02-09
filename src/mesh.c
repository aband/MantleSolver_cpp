#include "../include/mesh.h"
#include <petsc.h>

PetscErrorCode CreateFullMesh(DM dmMesh, Vec *fullmesh, MeshParam * mp){

    PetscErrorCode    ierr;
    Vec               gcoords;
    Vec               fmesh, lmesh;
    DM                dmCoords;
    PetscInt          xs,ys,xm,ym,M,N;
    Point             **localmesh;
    PetscInt          stencilwidth;
    PetscFunctionBeginUser;

    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dmMesh, &xs, &ys, NULL, &xm, &ym, NULL);                                                CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmMesh, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Define rectangular mesh
    ierr = DMGetLocalVector(dmMesh, &lmesh);                         CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dmMesh, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr); 
    ierr = DMGlobalToLocalEnd(dmMesh, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dmMesh, lmesh, &localmesh);               CHKERRQ(ierr);

    // Unpack parameters
    double L = mp->L;
    double H = mp->H;
    double xstart = mp->xstart;
    double ystart = mp->ystart;

    double hx = L/(double)M;
    double hy = H/(double)N;

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+xm; i++){
        localmesh[j][i].p[0] = xstart+i*hx;
        localmesh[j][i].p[1] = ystart+j*hy;
    }}

    ierr = DMDAVecRestoreArray(dmMesh, lmesh, &localmesh);       CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmMesh, &lmesh);                 CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
