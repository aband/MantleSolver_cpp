#ifndef TRANS_PARAM_H_
#define TRANS_PARAM_H_

#include "util.h"


// This is the only costumized function in simulation
double advfunc(const double& u, const vertex& vel, const vertex& unitnormal);

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work);

// Position function that returns where the cell is
std::string positin(const indice& gcell);

double InitCD(const valarray<double>& point, const vector<double>& param);

double InitHD(const valarray<double>& point, const vector<double>& param);

#endif
