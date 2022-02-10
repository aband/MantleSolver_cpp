#include "../include/output.h"

PetscErrorCode PrintFullMesh(DM dmMesh, Vec * fullmesh){

    PetscErrorCode ierr;
    Vec            fmesh, lmesh;
    PetscInt       xs,ys,xm,ym,M,N,stencilwidth;
    Point          **localmesh;
    PetscFunctionBeginUser;

    fmesh = *fullmesh; 

    ierr = DMDAGetCorners(dmMesh, &xs, &ys, NULL, &xm, &ym, NULL);                                                      CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmMesh, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Check Mesh definition
    DMGetLocalVector(dmMesh, &lmesh); 
    DMGlobalToLocalBegin(dmMesh,fmesh,INSERT_VALUES,lmesh);
    DMGlobalToLocalEnd(dmMesh,fmesh,INSERT_VALUES,lmesh);
    DMDAVecGetArrayRead(dmMesh,lmesh,&localmesh);

    for (int j=ys-stencilwidth; j<ym+ys+stencilwidth; j++){
    for (int i=xs-stencilwidth; i<xm+xs+stencilwidth; i++){
        printf("(%.2f,%.2f) ",localmesh[j][i].p[0],localmesh[j][i].p[1]);
    }printf("\n");}

    DMDAVecRestoreArrayRead(dmMesh,lmesh,&localmesh);
    DMRestoreLocalVector(dmMesh,&lmesh);

    PetscFunctionReturn(0);
}
