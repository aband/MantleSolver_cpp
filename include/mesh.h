#ifndef MESH_H_
#define MESH_H_

#include <time.h>
#include <petscdmda.h>

PetscErrorCode FullMesh();
PetscErrorCode TestControlMeshSecond();
PetscErrorCode TestControlMeshThird();
PetscErrorCode LogicRectMesh();
//PetscErrorCode CentroidGrid(); 
PetscErrorCode GetMeshInfo2D();

void DestroyMeshInfo2D();

#endif
