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


double basePoly(vertex& point, const vector<int>& param){
    return pow(point[0],param[0])*pow(point[1],param[1]);
}

// Evaluation of factorial

int factorial(int top, int bottom){
    assert(top>bottom || top==bottom);
    if (top==bottom){ return 1;}
    else{ return top*factorial(top-1,bottom);}
}

int factorial(int top){
    assert(top>0 || top==0);
    if (top==0){ return 1;}
    else{ return top*factorial(top-1);}
}

/*
 *int factorial(int n){
 *
 *    assert(n>0 || n==0);
 *
 *    int * work = new int [n+1] ();
 *
 *    work[0] = 1;
 *    for (int i=1; i<=n; i++){
 *        work[i] = i*work[i-1];
 *    }
 *    return work[n];
 *
 *}
 *
 *int factorial(int n, int m){
 *    assert(n>m || n==m);
 *
 *    int * work = new int [n-m] ();
 *
 *    if (n==m){
 *        return factorial(n);
 *    } else {
 *        work[0] = m+1;
 *        for (int i=1; i<n-m; i++){
 *            work[i] = (m+1+i)*work[i-1];
 *        }
 *        return work[n-m-1];
 *    }
 *}
 *
 */
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

// Compute factorial coefficient for polynomial derivatives
void polynDerMulti(int der, int max, int * multiplier){
    assert(der<max || der==max);

    for (int i=0; i<max-der; i++){
        multiplier[i] = factorial(der+i,i);
    }

}

indice MPILocalToGlobal(indice local, const MeshInfo& mi){
    return local + mi.MPIlocalCellStart;
}

indice MPIGlobalToLocal(indice global, const MeshInfo& mi){
    return global - mi.MPIlocalCellStart;
}

// Indice convention functions
// Flatten indice into 1D array
// Flatten into global indices!!!
int FlatIndic(const MeshInfo& mi, int i, int j)  
              {return j*mi.MPIglobalCellSize[0]+i;};

int FlatIndic(const int M, int i, int j) {return j*M+i;};
int FlatIndic(const MeshInfo& mi, const indice& p) {return FlatIndic(mi,p[0],p[1]);}
int FlatIndic(const int M, const indice& p) {return FlatIndic(M,p[0],p[1]);};

// Reverse process of flatten indices
indice Bend(const MeshInfo& mi, int flat) 
            {return {flat%mi.MPIglobalCellSize[0], flat/mi.MPIglobalCellSize[0]};};

indice Bend(const int M, int flat) {return {flat%M, flat/M};}

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

    // Pre calculate cell area for future computation.
    // Repeat calculation of cell areas cost a lot of computation resources.
    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+xm; i++){

        vertexSet corner;
        //! Retrieve local cell indice (including ghost vertex)
        indice ghostlayerShift {ghostWidth, ghostWidth};
        indice global {i,j};
        indice fullLocal = global - mi.MPIlocalCellStart + ghostlayerShift;

        for (auto & fcorner : mi.faceCorner){
            corner.push_back(mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0], fullLocal+fcorner)]); 
        }

        mi.cellArea.insert(std::make_pair<int,double>
                           (FlatIndic(mi.MPIglobalCellSize[0],i,j),
                            NumIntegralFace(corner,{0,0},{0.0,0.0},1.0,constFunc))); 
    }}

}

//! Extract corners for the target cell
vertexSet extractCorners(const MeshInfo& mi, const indice& global){
    vertexSet corner;

    //! Retrieve local cell indice (including ghost vertex)
    indice ghostlayerShift {mi.vertexGhostLayerSize, mi.vertexGhostLayerSize};
    indice fullLocal = global - mi.MPIlocalCellStart + ghostlayerShift;

    //! Extract corners from mesh.
    for (auto & fcorner: mi.faceCorner){
        corner.push_back(mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0],fullLocal+fcorner)]);
    }

    return corner;
}
