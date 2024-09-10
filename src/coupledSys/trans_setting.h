#ifndef TRANS_SETTING_H_
#define TRANS_SETTING_H_

// Return values of diffusion functions
// and advection functions
double diffFunc();

double advFunc();

// Return values on the boundary 
double bndryValDiff();

double bndryValAdv();

// Flux value prescribed on the boundary.
double bndryFluxAdv();

double bndryFluxDiff();

// Boundary categary


#endif
