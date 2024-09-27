#ifndef TRANS_PARAM_H_
#define TRANS_PARAM_H_

#include <array>
#include <string>
#include "util.h"

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

std::string location(const MeshInfo& mi,
                     const indice& globalCell);

#ifdef COUPLED
// Initial values
double InitHD(const vertex& point,
              const vector<double>& param);

double InitCD(const vertex& point,
              const vector<double>& param);
#endif

#endif
