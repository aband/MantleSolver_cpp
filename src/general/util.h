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

#include <functional>

#include "lapacke.h"
#include "integral.h"
#include <assert.h>

#include <petsc.h>

#define M_PI 3.14159265358979323846

using vertex = valarray<double>;
using indice = valarray<int>;

using vertexSet = vector<valarray<double>>;
using indiceSet = vector<valarray<int>>;

template <typename T>
using tensor = valarray<T>;

template <typename T>
using tensorSet = vector<valarray<T>>;

template <typename T>
void delete_pointed_to(T const ptr){
    delete ptr;
}

//! Print vertex and indice
template <typename T>
void nodePrint(T node){
    cout << "( " << node[0] << " ," << node[1] << " ) " ;
}

// ===== Define derivative class =====
using derivative = unordered_map<int, double>;

template <typename T>
T find_max(T val1, T val2){

    if (val1 > val2) {return val1;} else {return val2;}

}

template <typename T>
T harmonic_mean(T val1, T val2){

//    return 1.0/ (1.0/val1 + 1.0/val2);
    return 2*val1*val2/(val1+val2);
}

/**
 * Passing arithmetic function to the template 
 * std::minus
 * std::plus
 * std::multiplies
 * std::divides
 */
template <typename Key, typename Value, typename BinaryOp>
void unordered_map_arithmetic(std::unordered_map<Key, Value>& map1, 
                              const std::unordered_map<Key, Value>& map2, 
                                    BinaryOp operation) {

    for (auto it = map2.begin(); it != map2.end(); it++){
        auto findIt = map1.find(it->first);
        if (findIt != map1.end()){
            map1.at(it->first) = operation(map1.at(it->first),it->second);
        } else {
            map1.insert(std::pair<Key, Value> (it->first, operation(0.0,it->second)));
        }
    }
}

template<typename Key, typename Value, typename BinaryOp>
void unordered_map_arithmetic(std::unordered_map<Key, Value>& map,
                              const Value& scale,
                              BinaryOp operation){

    std::unordered_map<Key, Value> result;

    for (auto& pair : map){
        pair.second = operation(pair.second, scale);
    }
}

// A + xB style unordered_map arithmetic
template <typename Key, typename Value, typename BinaryOp1, typename BinaryOp2>
void unordered_map_arithmetic(std::unordered_map<Key, Value>& map1, 
                              const std::unordered_map<Key, Value>& map2, 
                                    BinaryOp1 operation1,
                              const Value& scale,
                                    BinaryOp2 operation2) {

    for (auto it = map2.begin(); it != map2.end(); it++){
        auto findIt = map1.find(it->first);
        if (findIt != map1.end()){
            map1.at(it->first) = operation1(map1.at(it->first),operation2(it->second,scale));
        } else {
            map1.insert(std::pair<Key, Value> (it->first, operation2(it->second,scale)));
        }
    }
}

template <typename Key, typename Value>
void unordered_map_print(const std::unordered_map<Key, Value>& map){

    for (auto it = map.begin(); it!=map.end(); it++){
        cout << std::setw(3) <<std::setprecision(3) << "( " << it->first <<  " " << it->second << " )   " ;
    } 
    cout << endl;
}

/**
 * Containing essential information about mesh.
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

    // Assign shift of ghost layer
    indice ghostShiftVertex;

    // Containing all the local mesh vertex points here  
    vertexSet lmesh; 

    // Corner index within a single element
    const indiceSet edgeCorner   {{0}, {1}};
    const indiceSet faceCorner   {{0,0},{1,0},{1,1},{0,1}};
    const indiceSet volumeCorner {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                                  {0,0,1},{1,0,1},{1,1,1},{0,1,1}};

    const indiceSet faceNormal {{0,-1},{1,0},{0,1},{-1,0}};

    // Scalar values distributed to mpi processors
    double** localVals;

    // Areas of cells distributed to mpi processors
    unordered_map<int,double> cellArea;

    // Coupled system only
    double** localCD;
    double** localHD;

    std::unordered_map<std::string , double**> localValsMap;

    // Physical length and height
    double L = 0.0;
    double H = 0.0;

} MeshInfo;

// Define function type of location functions
typedef bool (*LocFunc) (const indice& globalCell, 
                         const MeshInfo& mi);

void AssignValuesMeshInfo(MeshInfo& mi, DM dmv, DM dmu);

void printMeshInfo(MeshInfo& mi);

/**
 * Generic auxiliary functions.
 * Function return constant value.
 */
double constFunc(const valarray<double>& point,const vector<double>& param);

double constFunc(const valarray<double>& point,const vector<double>& param, double c);

double constFunc();

double constFunc(double c);

/**
 * Funcstions calculate factorials.
 */
int factorial(int top, int bottom);

int factorial(int top);

double basePoly(const vertex& point, const vector<int>& param);

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

//! Extract corners for the target cell
vertexSet extractCorners(const MeshInfo& mi, const indice& global);

/**!
 * Edges are arranged in the order that horizontal edges are numbered first,
 * while vertical edges are the next.
 * Three functions are defined as following representing the relationship 
 * between global cell index and global edge index.
 */
std::array<double,4> extractEdgeIndex(const MeshInfo& mi, const indice& globalCell);

std::array<indice,2> extractEdgeNbr(const MeshInfo& mi, const int& globalEdge);

std::array<vertex,2> extractEdge(const MeshInfo& mi, const int& globalEdge);

int edgeIndexLocalToGlobal(const MeshInfo& mi, const int& localEdge);
int edgeIndexGlobalToLocal(const MeshInfo& mi, const int& globalEdge);

//! Compute normalized tangent vector.
const vertex unitTangent(const vertexSet& edge, 
                         const double& len);

const vertex getUnitNormal(const std::array<vertex,2> edge,
                           const double& len);

const double getEdgeLength(const std::array<vertex, 2> edge);

// Some functions about cell and edge relationships

bool OutBndryCell(const MeshInfo& mi, 
                  const indice& gcell);

indice PickCellInside(const MeshInfo& mi,
                      const indice& gCellIn,
                      const indice& gCellOut);

template <typename T>
struct edgeEnds{
    T start;
    T end;
};

int extractVertEdgeInfo(const MeshInfo& mi, 
                        const indice& local,
                        const indice& ghostShift,
                        indice& gCellOut,
                        indice& gCellIn,
                        edgeEnds<vertex>& edgeEndsVertex,
                        edgeEnds<indice>& edgeEndsIndice);

int extractHoriEdgeInfo(const MeshInfo& mi,
                        const indice& local,
                        const indice& ghostShift,
                        indice& gCellOut,
                        indice& gCellIn,
                        edgeEnds<vertex>& edgeEndsVertex,
                        edgeEnds<indice>& edgeEndsIndice);

typedef int (*extractEdgeInfoFunc) (const MeshInfo& mi,
                                    const indice& local,
                                    const indice& ghostShift,
                                    indice& gCellOut,
                                    indice& gCellIn,
                                    edgeEnds<vertex>& edgeEndsVertex,
                                    edgeEnds<indice>& edgeEndsIndice);

#endif
