#ifndef TRANSFUNC_H_
#define TRANSFUNC_H_

#include "util.h"

// This is the only costumized function in simulation
// Transport functions
double advfunc(const double& u, const vertex& vel, const vertex& unitnormal);

double dfdu(const double& u);

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work);

#endif
