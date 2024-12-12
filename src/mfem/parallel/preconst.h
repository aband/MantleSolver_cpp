#ifndef PRECONST_H_
#define PRECONST_H_

#include <petsc.h>
#include <Hdivmixed.h>
#include <brmixed.h>
#include <util.h>
#include "passemble.h"

int SolScatAll(Vec *sol, Vec * g, 
               Vec * destSol, Vec * destg);

std::vector<double> GetFullSol(Vec * destSol, Vec * destg);

typedef struct{

    Vec vel_stokes, vel_darcy, g_stokes, g_darcy;

} ScatterResult;

// Create scatter map for all processor
// Repeated dof is fine here
template <typename T>
inline int scatterMap(T& funcSp, int m, int n, int M, int N,
                      PetscInt * lx          , PetscInt * ly,
                      PetscInt * scatMapInter, PetscInt * scatMapBndry,
                      std::unordered_map<int, int>& scatMapInterIndex,
                      std::unordered_map<int, int>& scatMapBndryIndex){


    // Create from index set
    int size = 0;

    if (m*n == 1) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "This is a sequential program \n"));
    } else {
        int yend   = 0;
        int xend   = 0;

        if (m == 1 ){
            // No partition in x direction
           for (int j=0; j<n-1; j++){
                yend += ly[j]; 
                // Loop through x direction
                for (int i=0; i<N; i++){
                    indice global {i, yend-1};
                    // face 3 is used by with another processor
                }
           } 
        }else if (n == 1){

        } else {


        } 

    }

    return 0;
}

template <typename T>
std::vector<double> ExtractWeightsParallel(PetscInt * lx, PetscInt * ly){

    // Scatter shared dofs ==========================
    VecScatter scatter; 

    // ==============================================

    std::vector<double> work;

    return work;
}

template <typename T>
int CGNSPrepareParallel(Vec * sol, Vec * g, 
                        const int *refmap, 
                        const MeshInfo& mi,
                        double * ux, double * uy,
                        T& funcSp, basis& basis_){

    CGNSPrepareParallel(sol, g, refmap, mi, ux, uy, funcSp, basis_, {0});

    return 1;
}

template <typename T>
int CGNSPrepareParallel(Vec * sol, Vec * g, 
                        const int *refmap, 
                        const MeshInfo& mi,
                        double * ux, double * uy,
                        T& funcSp, basis& basis_,
                        const std::vector<double>& parameter){

    int istart = mi.MPIlocalCellStart[0];
    int jstart = mi.MPIlocalCellStart[1];

    PetscScalar *valuesSol;
    PetscScalar *valuesg;

    VecGetArray(*sol, &valuesSol);
    VecGetArray(*g, &valuesg);

    for (int j=jstart; j<jstart + mi.MPIlocalCellSize[1]; j++){
    for (int i=istart; i<istart + mi.MPIlocalCellSize[0]; i++){

        indice global {i,j};

        vertex refcenter {0.0,0.0};

        basis_.GetCorners(mi, global);
        vertex center = GaussMapPointsFace(refcenter, basis_.corners());

        vertex val {0.0,0.0};

        std::vector<int> locdof = funcSp.LocalGlobalMap(mi, global);

        std::vector<vertex> basisVal = funcSp.EvaluateAll(basis_, center);

        indice local = global - mi.MPIlocalCellStart;
        int localflat = FlatIndic(mi.MPIlocalCellSize[0], local);

        for (int k=0; k<locdof.size(); k++){
//            if(funcSp.onBndry(mi, locdof.at(k))){
            if (bMarker(mi,funcSp.GlobalToLocalMapBndry(mi,locdof.at(k)),funcSp.name, parameter)==dirichlet){
                val += valuesg[refmap[locdof.at(k)]] * basisVal.at(k);
            } else {
                val += valuesSol[refmap[locdof.at(k)]] * basisVal.at(k);
            }
        }

        ux[localflat] = val[0];
        uy[localflat] = val[1];
    }}

    VecRestoreArray(*sol, &valuesSol);
    VecRestoreArray(*g, &valuesg);

    return 1;
}

// Compute velocity for a given cell with given local positions
template <typename T>
vector<vertex> ExtractVelocity(Vec * sol, Vec * g,
                               int *refmap,
                               const MeshInfo& mi,
                               vector<vertex> points,
                               const indice& gCell,
                               T& funcSp,
                               basis& mybasis){

   return ExtractVelocity(sol, g, refmap, mi, points, gCell, funcSp, mybasis, {0}); 
}

template <typename T>
vector<vertex> ExtractVelocity(Vec * sol, Vec * g,
                               int *refmap,
                               const MeshInfo& mi,
                               vector<vertex> points,
                               const indice& gCell,
                               T& funcSp,
                               basis& mybasis,
                               const std::vector<double>& parameter){

    std::vector<vertex> work;
    work.resize(points.size());

    PetscScalar *valuesSol;
    PetscScalar *valuesg;

    VecGetArray(*sol, &valuesSol);
    VecGetArray(*g, &valuesg);

    mybasis.GetCorners(mi, gCell);

    // !Get global indiex of the local dofs in specific and correct order
    const std::vector<int> elemDofs = funcSp.LocalGlobalMap(mi, gCell);

    for (int g=0; g<points.size(); g++){

        // Initialize interpolated value
        work.at(g) = {0.0,0.0};

        std::vector<vertex> basisVal = funcSp.EvaluateAll(mybasis, points.at(g));

        // Reconstruction of value with element basis
        for (int k=0; k<elemDofs.size(); k++){

            if (bMarker(mi,funcSp.GlobalToLocalMapBndry(mi,elemDofs.at(k)),funcSp.name, parameter) == dirichlet){
                work.at(g) += valuesg[refmap[elemDofs.at(k)]] * basisVal.at(k);
            } else {
                work.at(g) += valuesSol[refmap[elemDofs.at(k)]] * basisVal.at(k);
            }
        }
    }

    VecRestoreArray(*sol, &valuesSol);
    VecRestoreArray(*g, &valuesg);

    return work;
}

// Extract velocity on a given gauss points set
template <typename T>
int ExtractVelocityAll(unordered_map<int, vector<vertex>>& velocityAll,
                       const unordered_map<int, vector<vertex>>& edgeGaussPointsAll,
                       const MeshInfo& mi,
                       const int* refmap, Vec * sol, Vec * g, T& funcSp, basis& mybasis){

    indice local, gCellOut, gCellIn, gCell_inside;

    edgeEnds<vertex> edgeEndsVertex;
    edgeEnds<indice> edgeEndsIndice;

    vector<vertex> velocity_gaussp;

    for (const auto& [key, value] : edgeGaussPointsAll){

        if (key < mi.MPIlocalHoriEdgeSize){
            // Horizontal edge
            local = Bend(mi.MPIlocalCellSize[0], key);
            extractVertEdgeInfo(mi, local, mi.ghostShiftVertex, gCellOut, gCellIn, edgeEndsVertex, edgeEndsIndice);
        } else {
            // Vertical edge
            local = Bend(mi.MPIlocalVertexSize[0], key); 
            extractVertEdgeInfo(mi, local, mi.ghostShiftVertex, gCellOut, gCellIn, edgeEndsVertex, edgeEndsIndice);
        }

        gCell_inside = PickCellInside(mi, gCellIn, gCellOut); 
        velocity_gaussp = ExtractVelocity(sol, g, refmap, mi, value, gCell_inside,funcSp, mybasis);

        edgeGaussPointsAll.insert(
        make_pair<int, vector<vertex>>(key, velocity_gaussp));
 
    }

    return 1;
}

#endif
