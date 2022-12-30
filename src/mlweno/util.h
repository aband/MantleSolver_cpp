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

#include "lapacke.h"
#include "integral.h"
#include <assert.h>

using vertex = valarray<double>;
using indice = valarray<int>;

using vertexSet = vector<valarray<double>>;
using indiceSet = vector<valarray<double>>;

/*
 *Containing essential information about mesh.
 */
typedef struct {

    int dim;                 // Total number of dimensions

    vector<int> MPIlocalSize;         // Local chunck size without ghost layer
    indice MPIlocalStart;             // The starting vertex index for MPI local part

    vector<int> MPIglobalSize;        // global chunck size without ghost layer
    vector<int> cellGhostLayerSize;   // size of ghost layer of cell
    vector<int> vertexGhostLayerSize; // size of ghost layer of node

    // Horizontal edges first than vertical edges
    int MPIglobalHoriEdgeSize;
    int MPIlocalHoriEdgeSize;

    // Starting at 0 (has not stacked on horizontal edges yet)
    int MPIglobalVertEdgeSize;
    int MPIlocalVertEdgeSize;

    // Containing all the local mesh vertex points here  
    vertexSet lmesh; 
 
    // double** localVals;

    // Corner index within a single element
    const indiceSet edgeCorner   {{0}, {1}};
    const indiceSet FaceCorner   {{0,0},{1,0},{1,1},{0,1}};
    const indiceSet VolumeCorner {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                                     {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
} MeshInfo;

#endif
