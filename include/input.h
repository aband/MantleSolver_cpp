#ifndef INPUT_H_
#define INPUT_H_

#include <iostream>
#include <vector>
#include <valarray>
#include <petsc.h>

extern "C"{
#include "mesh.h"
}

using namespace std;

PetscErrorCode ReadMeshPortion(DM dm, Vec *fullmesh, vector< valarray<double> >& mesh);

#endif
