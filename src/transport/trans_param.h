#ifndef TRANS_PARAM_H_
#define TRANS_PARAM_H_

#include <array>

/**!
 * Transport boundary functions
 * No boundary condition function for flow conterparts for
 * the fact that flow boundary condition goes into assemble of linear system
 * and need to be compiled with mfem_lib
 */
double diffFunc();

std::array<double,2> advFunc(double u);

// Return values on the boundary 
std::array<double,2> bndryValDiff();

std::array<double,2> bndryValAdv();

// Flux value prescribed on the boundary.
double bndryFluxAdv();

double bndryFluxDiff();

#endif
