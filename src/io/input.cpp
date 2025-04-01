#include "input.h"

double constfun2(const valarray<double>& point,const vector<double>& param){
    return 1.0;
}

/*
 *This ReadMesh function transforms
 * an array in C to vector container in C++
 *holding mesh protion owned by current node.
 */

PetscErrorCode ReadMeshPortion(DM dm, Vec *fullmesh, vector< valarray<double> >& mesh){

    PetscErrorCode    ierr;
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

//This is a temperaroy function convert C style array to C++ style container
PetscErrorCode ReadSolutionLocal(DM dmu, Vec *globalu, vector< vector<double> >& sol){

    PetscErrorCode    ierr;
    Vec      gu = *globalu;
    Vec      lu;
    double   **localu; 
    PetscInt xs,ys,xm,ym,M,N,stencilwidth;
    PetscFunctionBeginUser;

    ierr = DMGetLocalVector(dmu, &lu);                      CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dmu, gu, INSERT_VALUES, lu);CHKERRQ(ierr); 
    ierr = DMGlobalToLocalEnd(dmu, gu, INSERT_VALUES, lu);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dmu, lu, &localu);               CHKERRQ(ierr);

    for (int j=ys-stencilwidth; j<ys+ym+stencilwidth; j++){
    for (int i=xs-stencilwidth; i<xs+xm+stencilwidth; i++){
         // This is specified for 2 dimension
         vector<double> target(localu[j][i],1);
         sol.push_back(target);
    }}

    ierr = DMDAVecRestoreArray(dmu, lu, &localu); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmu, &lu);           CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode SimpleInitialValue(DM dm, DM dmu, Vec *fullmesh, Vec *globalu, 
                                  double (*func)(const valarray<double>& point, const vector<double>& param)){

    PetscErrorCode ierr;
    Vec            gu,fmesh,lmesh;
    Point          **localmesh;
    double         **localu;
    PetscInt xs,ys,xm,ym,M,N,stencilwidth;
    PetscFunctionBeginUser;

    gu    = *globalu;
    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                       CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    ierr = DMDAVecGetArray(dmu,gu,&localu);CHKERRQ(ierr); 

    // Get local mesh
    ierr = DMGetLocalVector(dm, &lmesh);                         CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dm, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dm, lmesh, &localmesh);           CHKERRQ(ierr);

    Point   p0,p1,p2,p3;
    double  result;

    for (int j=ys; j<ym+ys; j++){
    for (int i=xs; i<xm+xs; i++){
        vector<valarray<double>> cell;

        p0 = localmesh[j][i+1];
        p1 = localmesh[j+1][i+1];
        p2 = localmesh[j+1][i];
        p3 = localmesh[j][i];

        valarray<double> point0(p0.p,2);
        cell.push_back(point0);
        valarray<double> point1(p1.p,2);
        cell.push_back(point1);
        valarray<double> point2(p2.p,2);
        cell.push_back(point2);
        valarray<double> point3(p3.p,2);
        cell.push_back(point3);

        vector<double> Empty;

        valarray<double> center = (cell[0]+cell[1]+cell[2]+cell[3])/4.0;

        double area = NumIntegralFace(cell,Empty,{0.0,0.0},1.0,constfun2);
        double result = NumIntegralFace(cell,{pow(area,0.5)},{0.0,0.0},1.0,func); 

        localu[j][i] = result/area;
    }}

    // Missing boundary condition
    // =====

    /* Restore local mesh */
    ierr = DMDAVecRestoreArrayRead(dm, lmesh, &localmesh);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &lmesh);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(dmu, gu, &localu);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode SimpleInitialValue(DM dm, DM dmu, Vec *fullmesh, Vec *globalu, const vector<double>& param,
                                  double (*func)(const valarray<double>& point, const vector<double>& param)){

    PetscErrorCode ierr;
    Vec            gu,fmesh,lmesh;
    Point          **localmesh;
    double         **localu;
    PetscInt xs,ys,xm,ym,M,N,stencilwidth;
    PetscFunctionBeginUser;

    gu    = *globalu;
    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                       CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    ierr = DMDAVecGetArray(dmu,gu,&localu);CHKERRQ(ierr); 

    // Get local mesh
    ierr = DMGetLocalVector(dm, &lmesh);                         CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dm, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dm, lmesh, &localmesh);           CHKERRQ(ierr);

    Point   p0,p1,p2,p3;
    double  result;

    for (int j=ys; j<ym+ys; j++){
    for (int i=xs; i<xm+xs; i++){
        vector<valarray<double>> cell;

        p0 = localmesh[j][i+1];
        p1 = localmesh[j+1][i+1];
        p2 = localmesh[j+1][i];
        p3 = localmesh[j][i];

        valarray<double> point0(p0.p,2);
        cell.push_back(point0);
        valarray<double> point1(p1.p,2);
        cell.push_back(point1);
        valarray<double> point2(p2.p,2);
        cell.push_back(point2);
        valarray<double> point3(p3.p,2);
        cell.push_back(point3);

        vector<double> Empty;

        valarray<double> center = (cell[0]+cell[1]+cell[2]+cell[3])/4.0;

        double area = NumIntegralFace(cell,Empty,{0.0,0.0},1.0,constfun2);
        double result = NumIntegralFace(cell,param,{0.0,0.0},1.0,func); 

        localu[j][i] = result/area;
    }}

    // Missing boundary condition
    // =====

    /* Restore local mesh */
    ierr = DMDAVecRestoreArrayRead(dm, lmesh, &localmesh);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &lmesh);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(dmu, gu, &localu);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode SimpleInitialValue(DM dm, Vec *globalu, const std::vector<double>& data){

    double **localu;
    Vec gu = *globalu;

    PetscInt xs,ys,xm,ym,M,N,stencilwidth;

    PetscCall(DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL)); 
    PetscCall(DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));

    PetscCall(DMDAVecGetArray(dm,gu,&localu));

    for (int j=ys; j<ym+ys; j++){
    for (int i=xs; i<xm+xs; i++){

        localu[j][i] = data.at(j*M+i);

    }}

    PetscCall(DMDAVecRestoreArray(dm, gu, &localu));

    return PETSC_SUCCESS;
}

PetscErrorCode ObliqueBurgers(DM dm, DM dmu, Vec *fullmesh, Vec *globalu, 
                              double (*func)(valarray<double>& point, const vector<double>& param)){

    PetscErrorCode ierr;
    Vec            gu,fmesh,lmesh;
    Point          **localmesh;
    double         **localu;
    PetscInt xs,ys,xm,ym,M,N,stencilwidth;
    PetscFunctionBeginUser;

    gu    = *globalu;
    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                       CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    ierr = DMDAVecGetArray(dmu,gu,&localu);CHKERRQ(ierr); 

    // Get local mesh
    ierr = DMGetLocalVector(dm, &lmesh);                         CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dm, fmesh, INSERT_VALUES, lmesh);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm, fmesh, INSERT_VALUES, lmesh);  CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dm, lmesh, &localmesh);           CHKERRQ(ierr);

    Point   p0,p1,p2,p3;
    double  result;

    for (int j=ys; j<ym+ys; j++){
    for (int i=xs; i<xm+xs; i++){
        vector<valarray<double>> cell;

        p0 = localmesh[j][i+1];
        p1 = localmesh[j+1][i+1];
        p2 = localmesh[j+1][i];
        p3 = localmesh[j][i];

        valarray<double> point0(p0.p,2);
        cell.push_back(point0);
        valarray<double> point1(p1.p,2);
        cell.push_back(point1);
        valarray<double> point2(p2.p,2);
        cell.push_back(point2);
        valarray<double> point3(p3.p,2);
        cell.push_back(point3);

        vector<double> Empty;

        valarray<double> center = (cell[0]+cell[1]+cell[2]+cell[3])/4.0;

        localu[j][i] = func(center, Empty);
    }}

    // Missing boundary condition
    // =====

    /* Restore local mesh */
    ierr = DMDAVecRestoreArrayRead(dm, lmesh, &localmesh);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &lmesh);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(dmu, gu, &localu);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
