#ifndef FUNC_H_
#define FUNC_H_

#include "reconstruction.h"

//! Two dimension functions
//! Functions for different dimensions are realized via function overloading.
double funcX(vertex x, double u, double t);

double dfuncX(vertex x, double u, double t);

double funcY(vertex x, double u, double t);

double dfuncY(vertex x, double u, double t);

#endif
