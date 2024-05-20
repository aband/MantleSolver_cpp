#ifndef error_H_
#define error_H_

#include <petsc.h>
#include "Hdivmixed.h"
#include "brmixed.h"
#include "util.h"
#include "myFunc.h"
#include "bndry.h"

// Compbine boundary values and computed solution
// Passed the test
// Correct output guaranteed
std::vector<double> GetFullSol(Vec * u, const bndryVal& bndryvals, int dof);

// Extract correct weights
std::array<double,8> ExtractWeights(const std::vector<double>& fullsol, 
                                    const std::array<int, 8> ltgMap);

std::array<double,12> ExtractWeights(const std::vector<double>& fullsol, 
                                    const std::array<int, 12> ltgMap);

// Return error measured in energy norm or any arbitrary norm
double L2ErrorElem(const std::array<double,8>& coeff,
                   const indice& globalElemIndic,
                   std::array<double,3> (*func)(const vertex& point),
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   Hdivmixed& hdiv_);

double L2ErrorElem(const std::array<double,12>& coeff,
                   const indice& globalElemIndic,
                   std::array<double,3> (*func)(const vertex& point),
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   BRMixed& br_);


double L2ErrorElem(const std::array<double,8>& coeff,
                   const indice& globalElemIndic,
                   vertex (*func)(const vertex& point, PhysProperty * pp),
                   PhysProperty * pp,
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   Hdivmixed& hdiv_);

double L2ErrorElem(const std::array<double,12>& coeff,
                   const indice& globalElemIndic,
                   vertex (*func)(const vertex& point, PhysProperty * pp),
                   PhysProperty * pp,
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   BRMixed& br_);

double L2ErrorElem(const double& approxP, 
                   std::array<double, 3> (*func)(const vertex& point),
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   const double& area);

/**
 * containing x, y coordinates and coresponding velocity vectors.
 * [x,y;vx,vy;exactx,exacty]
 * Data output to check error
 */
std::array<vertex, 3> quiverPrepare(const std::array<double, 12>& weight,
                                    const indice& globalElemIndic,
                                    const vertex& local,
                                    basis& basis_,
                                    BRMixed& br_,
                                    PhysProperty * pp);

std::array<vertex, 3> quiverPrepare(const std::array<double, 8>& weight,
                                    const indice& globalElemIndic,
                                    const vertex& local,
                                    basis& basis_,
                                    Hdivmixed& hdiv_,
                                    PhysProperty * pp);

int quiverOutput(const MeshInfo& mi, const std::vector<double>& fullSol, int M, int N, 
                 basis& basis_, BRMixed& br, Hdivmixed& hdiv, PhysProperty * pp, int flag);

int quiverOutput(const MeshInfo& mi, 
                 const std::vector<double>& fullSolStokes, 
                 const std::vector<double>& fullSolDarcy,
                 int M, int N,
                 basis& basis_, BRMixed& br, Hdivmixed& hdiv, PhysProperty * pp);

int quiverOutput(const MeshInfo& mi, 
                 const std::vector<double>& fullSolStokes, 
                 int M, int N,
                 basis& basis_, BRMixed& br, PhysProperty * pp);

// Parallel 
// Mixed interior and boundary dof cell
template <typename T>
inline std::vector<double> ExtractWeightsParallel(Vec * sol, Vec * g,
                                                  T& funcSp, const MeshInfo& mi,
                                                  const int * refmap,
                                                  const indice& global){

    std::vector<int> locdof = funcSp.LocalGlobalMap(mi, global); 

    std::vector<double> work;
    work.resize(locdof.size());

    Vec destination;

    VecCreateSeq(PETSC_COMM_SELF, 1 ,&destination);

    VecScatter scatter;
  
    IS from , to;
    PetscInt id_to = 0;

    PetscScalar *values;

    for (int i=0; i<locdof.size(); i++){
        const int id_from = refmap[locdof.at(i)];

        ISCreateGeneral(PETSC_COMM_SELF, 1, &id_from, PETSC_COPY_VALUES, &from);
        ISCreateGeneral(PETSC_COMM_SELF, 1, &id_to, PETSC_COPY_VALUES, &to);

        if (funcSp.onBndry(mi, locdof.at(i))){

            VecScatterCreate(*g, from, destination, to, &scatter);
            VecScatterBegin(scatter,*g, destination, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scatter,*g, destination, INSERT_VALUES, SCATTER_FORWARD);

        } else {

            VecScatterCreate(*sol, from, destination, to, &scatter);
            VecScatterBegin(scatter,*sol, destination, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scatter,*sol, destination, INSERT_VALUES, SCATTER_FORWARD);

        }

        VecGetArray(destination, &values);

        work.at(i) = values[0];
    }

    ISDestroy(&from);
    ISDestroy(&to);

    VecScatterDestroy(&scatter);

    return work;
};

// Create two arrays holding flow velocity at cell centroid 
template <typename T>
inline int CGNSPrepare(Vec * sol, Vec * g, 
                       const int * refmap,
                       const MeshInfo& mi,
                       double * ux, double * uy,
                       basis& basis_,
                       T& funcSp){

    int istart = mi.MPIlocalCellStart[0];
    int jstart = mi.MPIlocalCellStart[1];

    for (int j=jstart; j<jstart + mi.MPIlocalCellSize[1]; j++){
    for (int i=istart; i<istart + mi.MPIlocalCellSize[0]; i++){

        indice global {i,j};

        vertex refcenter {0.0,0.0};

        basis_.GetCorners(mi, global);
        vertex center = GaussMapPointsFace(refcenter, basis_.corners());

        vertex val {0.0,0.0};

        // Extract weight
        std::vector<double> weights = 
        ExtractWeightsParallel(sol, g, funcSp, mi, refmap, global);

        std::vector<vertex> work = funcSp.EvaluateAll(basis_, center);

        for (int loc=0; loc<work.size(); loc++){
            val += weights.at(loc)*work.at(loc); 
        }
  
        indice local = global - mi.MPIlocalCellStart;
        int localflat = FlatIndic(mi.MPIlocalCellSize[0], local);

        ux[localflat] = val[0];
        uy[localflat] = val[1];
    }}

    return 0;
}

#endif
