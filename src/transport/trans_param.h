#ifndef TRANS_PARAM_H_
#define TRANS_PARAM_H_

#include "util.h"

// This is the only costumized function in simulation
// Transport functions
double advfunc(const double& u, const vertex& vel, const vertex& unitnormal);

double dfdu(const double& u);

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work);

// Diffusion functions
double diffunc(const double& u);

// Position function that returns where the cell is
std::string location(const MeshInfo& mi, const indice& gcell);

std::string position(const indice& gcell);

double InitCD(const valarray<double>& point, const vector<double>& param);

double InitHD(const valarray<double>& point, const vector<double>& param);

#endif
