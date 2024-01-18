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

// Return error measured in energy norm or any arbitrary norm
double L2ErrorElem(const std::array<double,8>& coeff,
                   const indice& globalElemIndic,
                   std::array<double,3> (*func)(const vertex& point),
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   Hdivmixed& hdiv_);

double L2ErrorElem(const double& approxP, 
                   std::array<double, 3> (*func)(const vertex& point),
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   const double& area);

#endif
