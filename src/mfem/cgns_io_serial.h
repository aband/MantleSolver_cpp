#ifndef CGNS_IO_SERIAL_H_
#define CGNS_IO_SERIAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <petsc.h>
#include "util.h"
#include "error.h"

#include "cgnslib.h"

typedef struct {
    double x;
    double y;
} vec2D;

// Write mesh and cell centered flow
PetscErrorCode CgnsOutSerial(const MeshInfo& mi, const std::vector<double>& fullSol, 
                             int M, int N, basis& basis_, BRMixed& br, Hdivmixed& hdiv, 
                             PhysProperty * pp, int flag);


#endif
