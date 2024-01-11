#ifndef error_H_
#define error_H_

#include <petsc.h>
#include "Hdivmixed.h"
#include "brmixed.h"
#include "util.h"
#include "myFunc.h"
#include "bndry.h"

// Compbine boundary values and computed solution
std::vector<double> GetFullSol(Vec * u, const bndryVal& bndryvals, int dof);

// Extract correct weights
std::array<double,8> ExtractWeights(const std::vector<double>& fullsol, 
                                    const Hdivmixed& hdiv_,
                                    const indice& globalElem,
                                    const MeshInfo& mi);

// Return error measured in energy norm or any arbitrary norm
double L2ErrorElemInterior(const vector<double>& coeff,
                           const indice& globalElemIndic,
                           double (*func)(const vertex& point),
                           const valarray<double>& gwf,
                           const vector<vertex>& gpf,
                           basis& basis_,
                           Hdivmixed& hdiv_);

double L2ErrorElemBndry(const vector<double>& coeff,
                        const indice& globalElemIndic,
                        double (*func));

// Compute L2 error of divergence of u
// Overload function with different gauss points
double L2ErrorElemInterior();

double L2ErrorElemBndry();

#endif
