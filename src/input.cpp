#include "../include/input.h"

/*
 *This ReadMesh function transforms
 * an array in C to vector container in C++
 *holding mesh protion owned by current node.
 */

PetscErrorCode ReadMeshPortion(DM dm, Vec *fullmesh, vector< valarray<double> >& mesh){

    PetscErrorCode  ierr;
    Vec      fmesh = *fullmesh;  
    Vec      lmesh; 
    Point    **localmesh;
    PetscInt xs,ys,xm,ym,M,N,stencilwidth;
    PetscFunctionBeginUser; 

    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                                CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Define rectangular mesh
    ierr = DMGetLocalVector(dm, &lmesh);                         CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dm, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr); 
    ierr = DMGlobalToLocalEnd(dm, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dm, lmesh, &localmesh);               CHKERRQ(ierr);

    for (int j=ys-stencilwidth; j<ys+ym+stencilwidth; j++){
    for (int i=xs-stencilwidth; i<xs+xm+stencilwidth; i++){
         // This is specified for 2 dimension
         valarray<double> point(localmesh[j][i].p,2);
         mesh.push_back(point);
    }}

    ierr = DMDAVecRestoreArray(dm, lmesh, &localmesh); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &lmesh);           CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
