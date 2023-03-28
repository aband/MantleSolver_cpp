#ifndef MESH_H_
#define MESH_H_

#include <time.h>
#include <petscdm.h>
#include <petscsnes.h>

#define PI 3.141592653589793238

typedef struct {
    double xstart, ystart, L, H;
} MeshParam;

typedef struct {
    double p[2];
} Point;

PetscErrorCode CreateFullMesh(DM dmMesh, Vec * fullmesh, MeshParam * mp);
PetscErrorCode TestControlMeshSecond();
PetscErrorCode TestControlMeshThird();
PetscErrorCode LogicRectMesh(DM dmMesh, Vec * fullmesh, MeshParam * mp);
//PetscErrorCode CentroidGrid(); 
PetscErrorCode GetMeshInfo2D();

PetscErrorCode RefineMesh(DM dm, Vec * fullmesh, MeshParam * mp); 

void DestroyMeshInfo2D();

#endif
