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

PetscErrorCode PreconditionedUzawa(linearSys * ls, double tol, int MaxIter, double tau);

PetscErrorCode InexactUzawa(linearSys * ls, double tol, int MaxIter, double tau1, double tau2);

PetscErrorCode CoupledSolver(linearSys * ls1, linearSys * ls2, linearSys * lsResult, Mat * K, 
                             double tol, int MaxIter, double tau1, double tau2);

PetscErrorCode CoupledSolver(linearSys * ls1, linearSys * ls2, linearSys * lsResult, Mat * K, 
                             double tol, int MaxIter, double tau1);
#endif
