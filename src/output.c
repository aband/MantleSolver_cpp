#include "../include/output.h"

PetscErrorCode PrintFullMesh(DM dmMesh, Vec * fullmesh){

    PetscErrorCode ierr;
    Vec            fmesh, lmesh;
    PetscInt       xs,ys,xm,ym,M,N,stencilwidth;
    Point          **localmesh;
    PetscFunctionBeginUser;

    fmesh = *fullmesh; 

    ierr = DMDAGetCorners(dmMesh, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmMesh, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Check Mesh definition
    ierr = DMGetLocalVector(dmMesh, &lmesh); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dmMesh,fmesh,INSERT_VALUES,lmesh); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dmMesh,fmesh,INSERT_VALUES,lmesh); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dmMesh,lmesh,&localmesh); CHKERRQ(ierr);

    for (int j=ys-stencilwidth; j<ym+ys+stencilwidth; j++){
    for (int i=xs-stencilwidth; i<xm+xs+stencilwidth; i++){
        printf("(%.2f,%.2f) ",localmesh[j][i].p[0],localmesh[j][i].p[1]);
    }printf("\n");}

    ierr = DMDAVecRestoreArrayRead(dmMesh,lmesh,&localmesh); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmMesh,&lmesh); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode DrawPressure(DM dmu, Vec * globalu){

    PetscErrorCode    ierr;
    Vec      fu, lu;
    PetscInt xs,ys,xm,ym,M,N,stencilwidth;
    double   **localu;
    PetscFunctionBeginUser;

    fu = *globalu;

    ierr = DMDAGetCorners(dmu, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmu, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    FILE *f = fopen("Initial.data","w");

    ierr = DMGetLocalVector(dmu, &lu); CHKERRQ(ierr); 
    ierr = DMGlobalToLocalBegin(dmu,fu,INSERT_VALUES,lu); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dmu,fu,INSERT_VALUES,lu); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dmu,lu,&localu); CHKERRQ(ierr);

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+xm; i++){
        fprintf(f,"%f ",localu[j][i]);
    }fprintf(f,"\n");}

    ierr = DMDAVecRestoreArrayRead(dmu,lu,&localu); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmu,&lu); CHKERRQ(ierr);

    fclose(f);

    PetscFunctionReturn(0);
}
