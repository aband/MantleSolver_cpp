#include "error.h"

double assistTrueSolution(const vertex& point,
                          const vector<double>& param){
    return TrueSolution(point, param.at(0), param);
}

void L1Error(Vec *globalU,
             Vec *globalError,
             Vec *fullmesh,
             DM  dmu, 
             DM  dmMesh,
             const double& time,
             const MeshInfo& mi){

    int xs, ys, xm, ym;
    Vec localu, globalu; 
    Vec localE, globalE;

    DMDAGetCorners(dmu, &xs, &ys, NULL, &xm, &ym, NULL); 

    globalu = *globalU;
    globalE = *globalError;

    // Compute true cell averaged solution 
    SimpleInitialValue(dmMesh, dmu, fullmesh, &globalE, {time}, assistTrueSolution); 

    double **localuArray;
    double **localEArray;

    DMDAVecGetArray(dmu,globalu,&localuArray);
    DMDAVecGetArray(dmu,globalE,&localEArray);

    for (int j=ys; j<ym+ys; j++){
    for (int i=xs; i<xm+xs; i++){
        localEArray[j][i] = abs(localEArray[j][i] - 
	                             localuArray[j][i]);
    }}

    DMDAVecRestoreArray(dmu, globalu,&localuArray);
    DMDAVecRestoreArray(dmu, globalE,&localEArray);

}
