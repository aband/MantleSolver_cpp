#include "output.h"

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

    FILE *f = fopen("Pressure.data","w");

    ierr = DMGetLocalVector(dmu, &lu); CHKERRQ(ierr); 
    ierr = DMGlobalToLocalBegin(dmu,fu,INSERT_VALUES,lu); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dmu,fu,INSERT_VALUES,lu); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dmu,lu,&localu); CHKERRQ(ierr);

    for (int j=ys+ym/4; j<ys+ym-ym/4; j++){
    for (int i=xs+ym/4; i<xs+xm-ym/4; i++){
        fprintf(f,"%f ",localu[j][i]);
    }fprintf(f,"\n");}

    ierr = DMDAVecRestoreArrayRead(dmu,lu,&localu); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmu,&lu); CHKERRQ(ierr);

    fclose(f);

    PetscFunctionReturn(0);
}

//! Output data in a plain fashion.
PetscErrorCode PlainOutput(DM dmu, Vec * globalu, char* filename){

    PetscErrorCode    ierr;
    Vec      fu, lu;
    PetscInt xs,ys,xm,ym,M,N,stencilwidth;
    double   **localu;
    PetscFunctionBeginUser;

    fu = *globalu;

    ierr = DMDAGetCorners(dmu, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmu, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    FILE *f = fopen(filename,"w");

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

//! Output mesh in the plain form
PetscErrorCode PlainMeshOutput(DM dmMesh, Vec * fullmesh){

    PetscErrorCode ierr;
    Vec            fmesh, lmesh;
    PetscInt       xs,ys,xm,ym,M,N,stencilwidth;
    Point          **localmesh;
    PetscFunctionBeginUser;

    fmesh = *fullmesh; 

    ierr = DMDAGetCorners(dmMesh, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmMesh, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    FILE *f1 = fopen("gridX.txt", "w");
    FILE *f2 = fopen("gridY.txt", "w"); 

    // Check Mesh definition
    ierr = DMGetLocalVector(dmMesh, &lmesh); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dmMesh,fmesh,INSERT_VALUES,lmesh); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dmMesh,fmesh,INSERT_VALUES,lmesh); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dmMesh,lmesh,&localmesh); CHKERRQ(ierr);

    for (int j=ys; j<ym+ys+1; j++){
    for (int i=xs; i<xm+xs+1; i++){
        fprintf(f1,"%f ", localmesh[j][i].p[0]);
        fprintf(f2,"%f ", localmesh[j][i].p[1]);
    }fprintf(f1,"\n");fprintf(f2,"\n");}

    ierr = DMDAVecRestoreArrayRead(dmMesh,lmesh,&localmesh); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmMesh,&lmesh); CHKERRQ(ierr);

    fclose(f1);
    fclose(f2);

    PetscFunctionReturn(0);
}

// Draw matrix pattern
PetscErrorCode DrawMat(Mat V, const char * myfile){

    PetscFunctionBeginUser;

    FILE *f = fopen(myfile,"w");

    if (f == NULL){
        printf("Error opening file !\n");
        exit(1);
    }

    int mm,nn;
    MatGetSize(V,&nn,&mm);

    for (int j=nn-1; j>-1; j--){
    for (int i=0; i<mm; i++){
        double a;
        MatGetValues(V,1,&j,1,&i,&a);
        fprintf(f,"%f ",a);
    }fprintf(f,"\n ");}

    fclose(f);

    PetscFunctionReturn(0);
}

PetscErrorCode WriteMat(Mat V, const char * myfile){

    // Write matrix out in correct order for matlab

    PetscFunctionBeginUser;

    FILE *f = fopen(myfile, "w");

    if (f == NULL){
        printf("Error opening file !\n");
        exit(1);
    }

    int mm,nn;
    MatGetSize(V,&nn,&mm);

    for (int j=0; j<nn; j++){
    for (int i=0; i<mm; i++){
        double a;
        MatGetValues(V,1,&j,1,&i,&a);
        fprintf(f,"%.16lf ",a);
    }fprintf(f,"\n ");}

    fclose(f);

    PetscFunctionReturn(0);
}

PetscErrorCode WriteVec(Vec g, const char * myfile){

    PetscFunctionBeginUser;

    FILE *f = fopen(myfile, "w");

    if (f == NULL){
        printf("Error opening file !\n");
        exit(1);
    }

    int m;
    VecGetSize(g,&m);

    for (int i=0; i<m; i++){
        double a;
        VecGetValues(g,1,&i,&a);

        fprintf(f,"%.16lf ",a);
    }

    fclose(f);

    PetscFunctionReturn(0);
}

PetscErrorCode MPIIO(){
    PetscFunctionBeginUser;

    MPI_Comm comm = PETSC_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int err;

    char* filename = "testFile";
    MPI_File fh;

//    MPI_INT size;
//    MPI_INT rank;

//    MPI_Comm_size(PETSC_COMM_WORLD,&size);  
//    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);   

//    MPI_Offset offset = (MPI_Offset)rank*10*sizeof(int);

    int cmode;
    cmode = MPI_MODE_CREATE;
    cmode |= MPI_MODE_RDWR;

    err = MPI_File_open(comm, filename, cmode, info, &fh);
    assert(err == MPI_SUCCESS);


    MPI_File_close(&fh);

    PetscFunctionReturn(0);
}
