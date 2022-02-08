#include "../include/mesh.h"

PetscErrorCode FullMesh(DM dmCell, DM dmMesh, Vec *fullmesh, 
                        double L,      double H, 
                        double xstart, double ystart)
{
    PetscErrorCode    ierr;
    Vec               gcoords;
    Vec               fmesh, lmesh;
    DM                dmCoords;
    PetscInt          xs,ys,xm,ym,M,N;
    Vector2D          **coords, **localmesh;
    PetscInt          stencilwidth;
    PetscFunctionBeginUser;

    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dmMesh, &xs, &ys, NULL, &xm, &ym, NULL);                                                CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmMesh, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    ierr = DMGetCoordinateDM(dmCell, &dmCoords);       CHKERRQ(ierr);
    ierr = DMGetCoordinates(dmCell, &gcoords);         CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dmCoords, gcoords, &coords);CHKERRQ(ierr);

    // Define rectangular mesh
    ierr = DMGetLocalVector(dmMesh, &lmesh);                         CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dmMesh, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr); 
    ierr = DMGlobalToLocalEnd(dmMesh, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dmMesh, lmesh, &localmesh);               CHKERRQ(ierr);

    double hx = L/(double)M;
    double hy = H/(double)N;

    for (int j=ys-stencilwidth; j<ys+ym+stencilwidth; j++){
    for (int i=xs-stencilwidth; i<xs+xm+stencilwidth; i++){
        localmesh[j][i].x = xstart+i*hx;
        localmesh[j][i].y = ystart+j*hy;
    }}

    ierr = DMDAVecRestoreArray(dmMesh, lmesh, &localmesh);       CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(dmMesh,lmesh, ADD_VALUES, fmesh);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(dmMesh,lmesh,ADD_VALUES, fmesh);   CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmMesh, &lmesh);                 CHKERRQ(ierr);

    ierr = DMDAVecGetArray(dmMesh, fmesh, &localmesh);           CHKERRQ(ierr);

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+xm; i++){
        localmesh[j][i].x = coords[j][i].x;
        localmesh[j][i].y = coords[j][i].y;
    }}

    ierr = DMDAVecRestoreArray(dmCoords, gcoords, &coords);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(dmMesh, fmesh, &localmesh); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode TestControlMeshSecond(DM dmCell, double L, double H)
{
    PetscErrorCode    ierr;
    DM                dmCoords;
    Vec               globalcoords;
    Vector2D          **coords; 
    PetscInt          xs, ys, xm, ym;
    PetscFunctionBeginUser;

    ierr = DMGetCoordinateDM(dmCell, &dmCoords);                    CHKERRQ(ierr);
    ierr = DMGetCoordinates(dmCell, &globalcoords);                 CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dmCoords, globalcoords, &coords);        CHKERRQ(ierr);
    ierr = DMDAGetCorners(dmCoords, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    PetscMPIInt  size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    int      M, N;
    ierr = DMDAGetInfo(dmCell, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (M%4!=0 || N%4!=0){
        printf("Controlled mesh is only applied to even number of cells\n");
        PetscFunctionReturn(1);
    } else if (size != 1){
        printf("Controlled mesh will not run on parallel\n");
        PetscFunctionReturn(1);
    }

    double hy = H/(double)N;

    for (int j=1; j<N+1; j=j+2){
    for (int i=0; i<M  ; i++){
        coords[j][i].y = coords[j][i].y +hy*0.25*pow(-1,i%2+1);
    }}

    ierr = DMDAVecRestoreArray(dmCoords, globalcoords, &coords);CHKERRQ(ierr);

    PetscFunctionReturn(0);  
}

PetscErrorCode TestControlMeshThird(DM dmCell, double L, double H)
{
    PetscErrorCode    ierr;
    DM                dmCoords;
    Vec               globalcoords;
    Vector2D          **coords; 
    PetscInt          xs, ys, xm, ym;
    PetscFunctionBeginUser;

    ierr = DMGetCoordinateDM(dmCell, &dmCoords);                    CHKERRQ(ierr);
    ierr = DMGetCoordinates(dmCell, &globalcoords);                 CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dmCoords, globalcoords, &coords);        CHKERRQ(ierr);
    ierr = DMDAGetCorners(dmCoords, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    PetscMPIInt  size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    int      M, N;
    ierr = DMDAGetInfo(dmCell, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (M%4!=0 || N%4!=0){
        printf("Controlled mesh is only applied to even number of cells\n");
        PetscFunctionReturn(1);
    } else if (size != 1){
        printf("Controlled mesh will not run on parallel\n");
        PetscFunctionReturn(1);
    }

    double hy = H/(double)N;
    double hx = L/(double)M;

    double add[4] = {0,1,2,1};

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int iflag=i%4;
        switch(iflag)
        {
            case 1: coords[j][i].x = coords[j][i].x + hx*0.25*add[j%4]; break;
            case 3: coords[j][i].x = coords[j][i].x - hx*0.25*add[j%4]; break;
        }
        int jflag=j%4;
        switch(jflag)
        {
            case 1: coords[j][i].y = coords[j][i].y + hy*0.25*add[i%4]; break;
            case 3: coords[j][i].y = coords[j][i].y - hy*0.25*add[i%4]; break;
        } 
    }}

    ierr = DMDAVecRestoreArray(dmCoords, globalcoords, &coords);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode LogicRectMesh(DM dmCell, double L, double H)
{
    PetscErrorCode    ierr;
    DM                dmCoords;
    Vec               globalcoords;
    Vector2D          **coords; 
    PetscInt          xs, ys, xm, ym;
    PetscFunctionBeginUser;

    ierr = DMGetCoordinateDM(dmCell, &dmCoords);CHKERRQ(ierr);
    ierr = DMGetCoordinates(dmCell, &globalcoords);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dmCoords, globalcoords, &coords);CHKERRQ(ierr);
    ierr = DMDAGetCorners(dmCoords, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    int      M, N;
    ierr = DMDAGetInfo(dmCell, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    double   hx = L/(double)M, hy = H/(double)N;

    PetscRandom    rndx,rndy;
    PetscScalar    value;

    // Get current time in seconds as random seed
    time_t seconds;
    time(&seconds);

    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rndx);CHKERRQ(ierr);
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rndy);CHKERRQ(ierr);

    ierr = PetscRandomSetInterval(rndx,-0.3*hx,0.3*hx);  CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(rndx,(unsigned long)seconds);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rndx);                CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(rndy,-0.3*hy,0.3*hy);  CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(rndy,(unsigned long)seconds);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rndy);                CHKERRQ(ierr);

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+xm; i++){
        if(i*j !=0){
            ierr = PetscRandomGetValue(rndy, &value);CHKERRQ(ierr);
            coords[j][i].y = coords[j][i].y + (double)(PetscRealPart(value));
            ierr = PetscRandomGetValue(rndx, &value);CHKERRQ(ierr);
            coords[j][i].x = coords[j][i].x + (double)(PetscRealPart(value));
        }        
    }}

    ierr = DMDAVecRestoreArray(dmCoords, globalcoords, &coords);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode GetMeshInfo2D(DM dmMesh, Vec *fullmesh, MeshInfo2D ***meshinfo)
{
    PetscErrorCode    ierr; 
    Vec               fmesh, lmesh;
    Vector2D          **localmesh;
    PetscInt          xs, ys, xm, ym, M, N, s;
    PetscFunctionBeginUser;

    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dmMesh, &xs, &ys, NULL, &xm, &ym, NULL);                                     CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmMesh, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &s, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Memory allocation for meshinfo 
    MeshInfo2D *data     = (MeshInfo2D *)malloc((xm+2*s-1)*(ym+2*s-1)*sizeof(MeshInfo2D));
               *meshinfo = (MeshInfo2D **)malloc((ym+2*s-1)*sizeof(MeshInfo2D *));

    for (int i=0; i<ym+2*s-1; i++){
        (*meshinfo)[i] = &(data[(xm+2*s-1)*i]);
    }

    ierr = DMGetLocalVector(dmMesh, &lmesh);                         CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dmMesh, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dmMesh, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dmMesh, lmesh, &localmesh);           CHKERRQ(ierr);

    Vector2D p0,p1,p2,p3;

    for (int j=ys-s; j<ys+ym+s-1; j++){
    for (int i=xs-s; i<xs+xm+s-1; i++){
        p0 = localmesh[j][i];
        p1 = localmesh[j][i+1];
        p2 = localmesh[j+1][i+1];
        p3 = localmesh[j+1][i];
        (*meshinfo)[j-ys+s][i-xs+s].center.x = (p0.x+p1.x+p2.x+p3.x)/4.0;
        (*meshinfo)[j-ys+s][i-xs+s].center.y = (p0.y+p1.y+p2.y+p3.y)/4.0;
        (*meshinfo)[j-ys+s][i-xs+s].area     = 0.5*(PetscAbsReal((p2.x-p3.x)*(p0.y-p3.y) - (p2.y-p3.y)*(p0.x-p3.x)) + PetscAbsReal((p2.x-p1.x)*(p0.y-p1.y) - (p2.y-p1.y)*(p0.x-p1.x)) );
    }}

    ierr = DMDAVecRestoreArrayRead(dmMesh, lmesh, &localmesh);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmMesh, &lmesh);              CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

void DestoryMeshInfo2D(MeshInfo2D **meshinfo)
{
    int total = sizeof(meshinfo);
    int first = sizeof(meshinfo[0]);
    int NN = total/first;
    for (int i=0; i<NN; i++){
        free(meshinfo[i]);
    }
    free(meshinfo);
}
