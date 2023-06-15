#ifndef FUNC_H_
#define FUNC_H_

#include "reconstruction.h"

//! Two dimension functions
/**
 * Advection.
 */
double funcX(vertex x, double u, double t);

double dfuncX(vertex x, double u, double t);

double funcY(vertex x, double u, double t);

double dfuncY(vertex x, double u, double t);

/**
 * Diffusion.
 */
double diffFunc(double u);

double dDiffFunc(double u);

/**
 * Initial distribution.
 */
double InitialDistribution(vertex& point,
                           const vector<double>& param);

#endif
