#include "solve.h"

PetscErrorCode PreconditionedUzawa(linearSys * ls){

    /*
     * Using (preconditioned) CG for 
     * Using zero vectors as initial guesses.
     * x and y are zero vectors.
     */
    KSP ksp;
    PC  pc;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, ls->A, ls->A));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE)); // zero initial guess
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCBJACOBI));

    Mat BT;
    PetscCall(MatCreateTranspose(ls->B,&BT));

    // Derive x1 and y1
    Vec tmp1, tmp2;
    PetscCall(MatMult(BT, ls->y, tmp1));
    PetscCall(MatMult(ls->A, ls->x, tmp2));

    PetscScalar alpha = 1.0;
    PetscCall(VecAXPY(tmp2,alpha,tmp1));
    alpha = -1.0;
    PetscCall(VecAYPX(tmp2,alpha,ls->f));

    KSPSolve(ksp,tmp2,tmp1); 

    PetscCall(VecAXPY(ls->x,1,tmp1));



}
