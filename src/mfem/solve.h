#ifndef SOLVE_H_
#define SOLVE_H_

#include <petsc.h>
#include <iostream>
#include "bndry.h"

/*
 * Solve linear systems arised from mixed finite element scheme for 
 * Stokes and Darcy problems.
 *
 * Uzawa method solving for saddle point problem.
 * [A  -B^T] [x] = [f]
 * [-B -C  ] [y]   [-g]
 *
 * Preconditioned Uzawa solves A^{-1} exactly.
 * Inexact Uzawa solves A^{-1} exactly.
 */

typedef struct{
    Mat A, B, C;
    Vec f, g, x, y;

} linearSys; 

// Iterative solvers
PetscErrorCode CreateLinearSys(linearSys * ls, ReducedSys * reducedsys);

PetscErrorCode SimpleUzawa(linearSys * ls, double tol, int MaxIter, double tau, int precondType);

PetscErrorCode InexactUzawa(linearSys * ls, double tol, int MaxIter, double tau1, double tau2, int flag);

PetscErrorCode ExactUzawa(linearSys * ls, double tol, int MaxIter);

PetscErrorCode CoupledSolver(linearSys * ls1, linearSys * ls2, linearSys * lsResult, Mat * K, 
                             double tol, int MaxIter, double tau1, double tau2, int flag);

PetscErrorCode CoupledSolver(linearSys * ls1, linearSys * ls2, linearSys * lsResult, Mat * K, 
                             double tol, int MaxIter, double tau1);

PetscErrorCode CoupledExactUzawa(linearSys * ls, double tau1, double tau2, 
                                 double tol, int MaxIter, int pType);

#endif
