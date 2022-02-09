#ifndef MESH_H_
#define MESH_H_

#include <time.h>
#include <petscdm.h>
#include <petscsnes.h>

PetscErrorCode CreateFullMesh();
PetscErrorCode TestControlMeshSecond();
PetscErrorCode TestControlMeshThird();
PetscErrorCode LogicRectMesh();
//PetscErrorCode CentroidGrid(); 
PetscErrorCode GetMeshInfo2D();

void DestroyMeshInfo2D();

typedef struct {
    double xstart, ystart, L, H;
} MeshParam;

typedef struct {
    double p[2];
} Point;

#endif
