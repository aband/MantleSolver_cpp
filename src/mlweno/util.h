#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <utility>
#include <set>
#include <unordered_set>
#include <array>
#include <valarray>
#include <algorithm>
#include <numeric>
#include <memory>
#include <cstdlib>
#include <type_traits>
#include <iomanip>

#include "lapacke.h"
#include "integral.h"
#include <assert.h>

#include <petsc.h>

using vertex = valarray<double>;
using indice = valarray<int>;

using vertexSet = vector<valarray<double>>;
using indiceSet = vector<valarray<double>>;

template <typename T>
void delete_pointed_to(T const ptr){
    delete ptr;
}

/*
 *Containing essential information about mesh.
 */
typedef struct {

    //! In the case, Cell and Vertex are maintained by same global size
    //! They will still be stored separately for clearification.

    int dim;                 //! Total number of dimensions

    vector<int> MPIlocalCellSize;         //! Local chunk size of cell without ghost layer
    vector<int> MPIlocalVertexSize;       //! Local chunk size of vertex without ghost layer

    indice MPIlocalCellStart;             //! The starting cell index for MPI local part
    indice MPIlocalVertexStart;           //! The starting vertex index for MPI local part

    vector<int> MPIglobalCellSize;        //! global chunk size of cell without ghost layer
    vector<int> MPIglobalVertexSize;      //! global chunk size of vertex without ghost layer

    int cellGhostLayerSize;   //! size of ghost layer of cell
    int vertexGhostLayerSize; //! size of ghost layer of node

    vector<int> MPIlocalCellSizeFull;     //! Local chunk size of cells including ghost layer
    vector<int> MPIlocalVertexSizeFull;   //! Local chunk size of vertex including ghost layer 

    // Horizontal edges first than vertical edges
    int MPIglobalHoriEdgeSize;
    int MPIlocalHoriEdgeSize;

    // Starting at 0 (has not stacked on horizontal edges yet)
    int MPIglobalVertEdgeSize;
    int MPIlocalVertEdgeSize;

    // Containing all the local mesh vertex points here  
    vertexSet lmesh; 
 
    double** localVals;

    // Corner index within a single element
    const indiceSet edgeCorner   {{0}, {1}};
    const indiceSet faceCorner   {{0,0},{1,0},{1,1},{0,1}};
    const indiceSet volumeCorner {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                                  {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
} MeshInfo;

void AssignValuesMeshInfo(MeshInfo& mi, DM dmv, DM dmu);

// Generic auxiliary functions
// Function return constant value
double constFunc(valarray<double>& point,const vector<double>& param);

double constFunc(valarray<double>& point,const vector<double>& param, double c);

double constFunc();

double constFunc(double c);

// Funcstions calculate factorials
int factorial(int top, int bottom);

int factorial(int top);

double basePoly(vertex& point, const vector<int>& param);

// Evaluation of polynomial using Horner's method
double polyEval(double x, double * coef, int degree);

// Compute factorial coefficient for polynomial derivatives
void polynDerMulti(int der, int max, int * multiplier);

// Local to global and global to local
// All indices are referenced to cell indice
indice MPILocalToGlobal(indice local, const MeshInfo& mi);
indice MPIGlobalToLocal(indice global, const MeshInfo& mi);

// Indice convention functions
// Flatten indice into 1D array
int FlatIndic(const MeshInfo& mi, int i, int j);

int FlatIndic(const int M, int i, int j);
int FlatIndic(const MeshInfo& mi, const indice& p);
int FlatIndic(const int M, const indice& p);

// Reverse process of flatten indices
indice Bend(const MeshInfo& mi, int flat);

indice Bend(const int M, int flat);

#endif
