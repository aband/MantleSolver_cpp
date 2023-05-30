#ifndef CGNS_IO_H_
#define CGNS_IO_H_

#include <stdio.h>
#include <stdlib.h>
#include <petsc.h>

#include "pcgnslib.h"
#include "mpi.h"

typedef struct {
  double   x;
  double   y;
} Vector2D;

//PetscErrorCode DMDACgnsOut2D(DM dmMesh, Vec * fullmesh, DM dmCell, Vec * Sol, char * filename);

PetscErrorCode CGNSMeshWrite(DM dm, Vec * fullmesh);

#endif
