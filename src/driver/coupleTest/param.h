#ifndef PARAM_H_
#define PARAM_H_

#include "util.h"
#include "eutectic.h"
#include "myFunc.h"

using namespace EUTECTIC;

// Initial distribution of c_bar and h_bar
double ComputePorosity(const vertex& point, PhysProperty * pp, phaseState * pPtr);

int PorosityOut(double xstart, double ystart, double L, double H, int seed, 
                PhysProperty * pp);

#endif
