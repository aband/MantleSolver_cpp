#include "util.h"

// Couple of constant functions
double constFunc(valarray<double>& point,const vector<double>& param){
    return 1.0;
}

double constFunc(valarray<double>& point,const vector<double>& param, double c){
    return c;
}

double constFunc(){
    return 1.0;
}

double constFunc(double c){
    return c;
}


// Evaluation of factorial
/*
 *int factorial(int top, int bottom){
 *    assert(top>bottom || top==bottom);
 *    if (top==bottom){ return 1;}
 *    else{ return top*factorial(top-1,bottom);}
 *}
 */

/*
 *int factorial(int top){
 *    assert(top>0 || top==0);
 *    if (top==0){ return 1;}
 *    else{ return top*factorial(top-1);}
 *}
 */

int factorial(int n){

    assert(n>0 || n==0);

    int * work = new int [n+1] ();

    work[0] = 1;
    for (int i=1; i<=n; i++){
        work[i] = i*work[i-1];
    }
    return work[n];

}

int factorial(int n, int m){
    assert(n>m || n==m);

    int * work = new int [n-m] ();

    work[0] = m+1;
    for (int i=1; i<n-m; i++){
        work[i] = (m+1+i)*work[i-1];
    }
    return work[n-m-1];
}

// Evalutaion of polynomial using Horner's method
double polyEval(double x, double * coef, int degree){

    if (abs(x) <= 1){

        double work = coef[degree];

        for (int r=degree-1; r>=0; r--){
            work = work*x + coef[r]; 
        } 

        return work;

    } else {

        double work = coef[0];

        for (int r=1; r<=degree; r++){
            work = work/x + coef[r];
        }

        return pow(x,degree)*work;
    }

}

// Assign values to MeshInfo object
void AssignValuesMeshInfo(MeshInfo& mi, DM dmv, DM dmu){

    PetscInt     dim, xs, ys, xm, ym, M, N;
    PetscInt     ghostWidth;

    PetscFunctionBeginUser;

    // Extract information of solution u
    DMDAGetCorners(dmu, &xs, &ys, NULL, &xm, &ym, NULL);
    DMDAGetInfo(dmu, &dim, &M, &N, NULL, NULL, NULL, NULL, NULL, &ghostWidth, NULL, NULL, NULL, NULL);

    // Assign values to meshInfo members based on the information above
    mi.MPIlocalCellSize.push_back(xm);
    mi.MPIlocalCellSize.push_back(ym);

    mi.MPIlocalCellStart = {xs,ys};

    mi.MPIglobalCellSize.push_back(M);
    mi.MPIglobalCellSize.push_back(N);
    mi.cellGhostLayerSize = ghostWidth; 

    mi.MPIlocalCellSizeFull.push_back(xm+2*ghostWidth);
    mi.MPIlocalCellSizeFull.push_back(ym+2*ghostWidth);

    // Extract information of vertex dm
    DMDAGetCorners(dmv, &xs, &ys, NULL, &xm, &ym, NULL);
    DMDAGetInfo(dmv, &dim, &M, &N, NULL, NULL, NULL, NULL, NULL, &ghostWidth, NULL, NULL, NULL, NULL);

    mi.MPIlocalVertexSize.push_back(xm);
    mi.MPIlocalVertexSize.push_back(ym);

    mi.MPIlocalVertexStart = {xs,ys};

    mi.MPIglobalVertexSize.push_back(M);
    mi.MPIglobalVertexSize.push_back(N);
    mi.vertexGhostLayerSize = ghostWidth; 

    mi.MPIlocalVertexSizeFull.push_back(xm+2*ghostWidth);
    mi.MPIlocalVertexSizeFull.push_back(ym+2*ghostWidth);

}
