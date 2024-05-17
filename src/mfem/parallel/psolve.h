#ifndef PSOLVE_H_
#define PSOLVE_H_

#include <petsc.h>
#include <iostream>
#include "passemble.h"

// parallel uzawa iteration

int CreateLinearSys(ReducedSys * redsys, const int& nelem);

int CreateCoupledSystem(ReducedSys * redsys1, ReducedSys * redsys2, 
                        ReducedSys * redsysRes, Mat * K);

int CoupledUzawa(ReducedSys * redsys, double tol, int MaxIter);

#endif
