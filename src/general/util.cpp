#include "util.h"

// Couple of constant functions
double constFunc(const valarray<double>& point,const vector<double>& param){
    return 1.0;
}

double constFunc(const valarray<double>& point,const vector<double>& param, double c){
    return c;
}

double constFunc(){
    return 1.0;
}

double constFunc(double c){
    return c;
}


double basePoly(const vertex& point, const vector<int>& param){
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

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Extract information of vertex dm
    DMDAGetCorners(dmv, &xs, &ys, NULL, &xm, &ym, NULL);
    DMDAGetInfo(dmv, &dim, &M, &N, NULL, NULL, NULL, NULL, NULL, &ghostWidth, NULL, NULL, NULL, NULL);

    mi.MPIlocalVertexSize.push_back(xm+1);
    mi.MPIlocalVertexSize.push_back(ym+1);

    mi.MPIlocalVertexStart = {xs,ys};

    mi.MPIglobalVertexSize.push_back(M+1);
    mi.MPIglobalVertexSize.push_back(N+1);
    mi.vertexGhostLayerSize = ghostWidth; 

    mi.MPIlocalVertexSizeFull.push_back(xm+2*ghostWidth);
    mi.MPIlocalVertexSizeFull.push_back(ym+2*ghostWidth);

    // Assign values to edge number
    mi.MPIglobalHoriEdgeSize = mi.MPIglobalCellSize[0]*mi.MPIglobalVertexSize[1];
    mi.MPIglobalVertEdgeSize = mi.MPIglobalCellSize[1]*mi.MPIglobalVertexSize[0];

    // Pre calculate cell area for future computation.
    // Repeat calculation of cell areas cost a lot of computation resources.
    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    for (int j=ys; j<ys+ym+1; j++){
    for (int i=xs; i<xs+xm+1; i++){


        if (j<N && i<M){
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
            }
    }}

}

void printMeshInfo(MeshInfo& mi){
    cout << "local size of Cells :" << mi.MPIlocalCellSize[0] << " " << mi.MPIlocalCellSize[1] << endl;

    cout << "local size of Vertexs :" << mi.MPIlocalVertexSize[0] << " " << mi.MPIlocalVertexSize[1] << endl;

    cout << "local start of Cells :" << mi.MPIlocalCellStart[0] << " " << mi.MPIlocalCellStart[1] << endl;

    cout << "local start of Vertexs :" << mi.MPIlocalVertexStart[0] << " " << mi.MPIlocalVertexStart[1] << endl;


    cout << "global size of cells :" << mi.MPIglobalCellSize[0]<< " " << mi.MPIglobalCellSize[1] << endl;

    cout << "global size of vertexs :" << mi.MPIglobalVertexSize[0]<< " " << mi.MPIglobalVertexSize[1] << endl;

    cout << "global size of vertical edges :" << mi.MPIglobalVertEdgeSize << endl;
    cout << "global size of horizontal edges :" << mi.MPIglobalHoriEdgeSize << endl;

    cout << endl;
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

// =============================================================================
// local edges are ordered as following:
//   __3___
//  |      |
// 0|      |2
//  |______|
//     1
// =============================================================================
std::array<double,4> extractEdgeIndex(const MeshInfo& mi, const indice& globalCell){
    std::array<double, 4> work;

    // Two horizontal edges counted first.
    work[1] = FlatIndic(mi, globalCell); 

    work[3] = FlatIndic(mi, {globalCell[0], globalCell[1]+1});

    work[0] = mi.MPIglobalHoriEdgeSize + FlatIndic(mi.MPIglobalCellSize[0]+1, globalCell);

    work[2] = mi.MPIglobalHoriEdgeSize + FlatIndic(mi.MPIglobalCellSize[0]+1, 
                                         {globalCell[0]+1,globalCell[1]});

    return work;
}

// Return two neighbours of this given edge index.
std::array<indice,2> extractEdgeNbr(const MeshInfo& mi, const int& globalEdge){
    std::array<indice,2> nBr;

    if (globalEdge < mi.MPIglobalHoriEdgeSize){
        // Horizontal edge
        nBr[0] = Bend(mi,globalEdge);
        nBr[1] = {nBr.at(0)[0], nBr.at(0)[1]-1};
    } else {
        // Vertical edge
        nBr[0] = Bend(mi.MPIglobalCellSize[0]+1,
                      globalEdge-mi.MPIglobalHoriEdgeSize);
        nBr[1] = {nBr.at(0)[0]-1, nBr.at(0)[1]};
    }

    return nBr;
}
/**
 * Extract begin and end points for the given global edge index.
 * Two points are ordered the same,
 * From left to right for horizontal edge
 * From bottom to top for vertical edge
 */
std::array<vertex,2> extractEdge(const MeshInfo& mi, const int& globalEdge){

    std::array<vertex,2> edge;

    indice ghostlayerShift {mi.vertexGhostLayerSize, mi.vertexGhostLayerSize};

    if (globalEdge < mi.MPIglobalHoriEdgeSize){
        indice leftVertexIndex = Bend(mi,globalEdge) - mi.MPIlocalCellStart + 
                                 ghostlayerShift;

        edge[0] = mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0], leftVertexIndex)];

        leftVertexIndex += {1,0};

        edge[1] = mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0], leftVertexIndex)];

    } else {
        indice bottomVertexIndex = Bend(mi.MPIglobalVertexSize[0], globalEdge) - mi.MPIlocalCellStart + ghostlayerShift;

        edge[0] = mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0],
                  bottomVertexIndex)];

        bottomVertexIndex += {0,1};

        edge[1] = mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0],
                  bottomVertexIndex)];
    }

    return edge;
}

// Return global edge index corresponding to local index
// ghost region not included.
int edgeIndexGlobalToLocal(const MeshInfo& mi,
                           const int& globalEdge){

    std::array<indice,2> gcells = extractEdgeNbr(mi, globalEdge);

    indice localcell = MPIGlobalToLocal(gcells[0],mi);
   
    int work = 0;

    if (globalEdge < mi.MPIglobalHoriEdgeSize){
        // Horizontal edge
        work = FlatIndic(mi.MPIlocalCellSize[0], localcell);
    } else {
        // Vertical edge
        work = FlatIndic(mi.MPIlocalCellSize[0]+1, localcell);
    }
  
    return work;
}

int edgeIndexLocalToGlobal(const MeshInfo& mi,
                           const int& localEdge){

    indice localCell;

    if (localEdge < mi.MPIlocalHoriEdgeSize){
        // It is a Horizontal edge
        localCell = Bend(mi.MPIlocalCellSize[0], localEdge);
        return FlatIndic(mi, localCell+mi.MPIlocalCellStart); 
    } else {
        // It is a Vertical edge
        localCell = Bend(mi.MPIlocalCellSize[0]+1, localEdge - mi.MPIlocalVertEdgeSize);
        return FlatIndic(mi, localCell+mi.MPIlocalCellStart) + mi.MPIglobalHoriEdgeSize;
    }
}

// =============================================================================

const vertex unitTangent(const vertexSet& edge, const double& len){
    return (edge.at(1) - edge.at(0))/len;
}

const vertex getUnitNormal(const std::array<vertex, 2> edge, 
                           const double& len){

    vertex work;
    work = edge[1] - edge[0];
    work = work.cshift(1);
    work[1] *= -1;
    work /= len;

    return work;
}

const double getEdgeLength(const std::array<vertex, 2> edge){

    vertex vec = edge[1]-edge[0];
    vec *= vec;
    return sqrt(vec.sum());
}
