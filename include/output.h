#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <stdio.h>
#include <stdlib.h>

#include <petscdmda.h>
#include "../include/mesh.h"

PetscErrorCode PrintFullMesh(DM dmMesh, Vec * fullmesh);

#endif
