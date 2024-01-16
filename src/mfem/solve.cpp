#include "solve.h"

PetscErrorCode PreconditionedUzawa(linearSys * ls, double tol, int MaxIter, double tau){

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

    Mat B;
    PetscCall(MatCreateTranspose(ls->B,&B));

    double r = 1.0;
    int    iter = 0;

    int cM, cN;

    VecGetSize(ls->f, &cM);
    VecGetSize(ls->g, &cN);

    Vec tmp1, tmp2, tmp3;

    PetscCall(VecCreate(PETSC_COMM_WORLD, &tmp1));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &tmp2));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &tmp3));

    PetscCall(VecSetSizes(tmp1,PETSC_DECIDE,cM));
    PetscCall(VecSetSizes(tmp2,PETSC_DECIDE,cM));
    PetscCall(VecSetSizes(tmp3,PETSC_DECIDE,cN));

    PetscCall(VecSetUp(tmp1));
    PetscCall(VecSetUp(tmp2));
    PetscCall(VecSetUp(tmp3));

    PetscCall(VecZeroEntries(tmp1));
    PetscCall(VecZeroEntries(tmp2));
    PetscCall(VecZeroEntries(tmp3));

    while(r > tol && iter < MaxIter){

        PetscCall(MatMult(ls->B, ls->y, tmp1));
        PetscCall(MatMult(ls->A, ls->x, tmp2));

        PetscScalar alpha = 1.0;
        PetscCall(VecAXPY(tmp2,alpha,tmp1));
        alpha = -1.0;
        PetscCall(VecAYPX(tmp2,alpha,ls->f));

        KSPSolve(ksp,tmp2,tmp1); 

        PetscCall(VecAXPY(ls->x,1,tmp1)); // x1

        PetscCall(MatMult(B,ls->x,tmp3));
        PetscCall(VecAXPY(tmp3, -1, ls->g));
        PetscCall(VecAXPY(ls->y,tau,tmp3));

        PetscReal val1, val2;
        PetscCall(VecNorm(tmp1,NORM_2,&val1));
        PetscCall(VecNorm(tmp3,NORM_2,&val2));
        r = val1 + val2; 

        iter++;
    }

    if (iter < MaxIter){
        printf("Uzawa converged successfully! r = %.3e, Used %d iterations. \n", r, iter);
        return PETSC_SUCCESS;
    } else {
        printf("Uzawa failed to converge! r = %.3e \n", r);
        return PETSC_ERR_CONV_FAILED;
    }
}
