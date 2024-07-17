#ifndef PARAM_H_
#define PARAM_H_

#include "util.h"

// Initial distribution of c_bar and h_bar
double NDcompbar(const vertex& point, 
                 const double& phi_f,
                 const double& c_bar,
                 PhysProperty * pp);

double NDhbar(const vertex& point);

#endif
