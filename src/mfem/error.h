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

#endif
