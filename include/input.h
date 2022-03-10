#ifndef INPUT_H_
#define INPUT_H_

#include <iostream>
#include <vector>
#include <valarray>
#include <petsc.h>
#include <numeric>

#include "integral.h"

extern "C"{
#include "mesh.h"
}

using namespace std;

PetscErrorCode ReadMeshPortion(DM dm, Vec *fullmesh, vector< valarray<double> >& mesh);
PetscErrorCode SimpleInitialValue(DM dm, DM dmu, Vec *fullmesh, Vec *globalu, 
                                  double (*func)(valarray<double>& point, const vector<double>& param)); 
PetscErrorCode ReadSolutionLocal(DM dmu, Vec *globalu, vector< vector<double> >& sol);


#endif
