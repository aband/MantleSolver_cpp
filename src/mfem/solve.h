#ifndef SOLVE_H_
#define SOLVE_H_

/*
 * Solve linear systems arised from mixed finite element scheme for 
 * Stokes and Darcy problems.
 *
 * Uzawa method solving for saddle point problem.
 * [A B^T] [x] = [f]
 * [B O  ] [y]   [g]
 *
 * Preconditioned Uzawa solves A^{-1} exactly.
 * Inexact Uzawa solves A^{-1} exactly.
 */

typedef struct{
    Mat A, B, C;
    Vec f, g, x, y;

} linearSys; 

PetscErrorCode PreconditionedUzawa(Mat * A, Mat * B, 
                                   Vec * f, Vec * g, 
                                   Vec * x, Vec * y){


}

PetscErrorCode InexectUzawa(Mat * A, Mat * B, Vec * f, Vec * g, Vec * x, Vec * y){


}

#endif
