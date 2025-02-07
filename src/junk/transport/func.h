#ifndef FUNC_H_
#define FUNC_H_

//#include "reconstMLWENO.h"
#include "reconstruction.h"

//! Two dimension functions
/**
 * Advection.
 */
double funcX(vertex x, double u, double t);

double dfuncX(vertex x, double u, double t);

double funcY(vertex x, double u, double t);

double dfuncY(vertex x, double u, double t);

double funcX(double u);

double dfuncX(double u);

double funcY(double u);

double dfuncY(double u);

/**
 * Diffusion.
 */
double diffFunc(double u);

double dDiffFunc(double u);

/**
 * Initial distribution.
 */
double InitialDistribution(const vertex& point,
                           const vector<double>& param);

/**
 * True solution.
 */

double TrueSolution(const vertex& point, 
                    const double& time,
                    const vector<double>& param); 

#endif
