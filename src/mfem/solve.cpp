#include "solve.h"
#include <iostream>

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
    PetscCall(KSPCGSetType(ksp, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE)); // zero initial guess
    PetscCall(KSPGetPC(ksp, &pc));
    //PetscCall(PCSetType(pc, PCBJACOBI));

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

    Vec test;
    PetscCall(VecCreate(PETSC_COMM_WORLD, &test));
    PetscCall(VecSetSizes(test,PETSC_DECIDE,cM));
    PetscCall(VecSetUp(test));
    PetscCall(VecZeroEntries(test));

    while(r > tol && iter < MaxIter){

        std::cout << "Iter : " << iter+1 << std::endl << std::endl;

        PetscCall(MatMult(ls->B, ls->y, tmp1));
        PetscCall(MatMult(ls->A, ls->x, tmp2));

        PetscCall(VecAXPY(tmp2,1.0,tmp1));

        PetscCall(VecAYPX(tmp2,-1.0,ls->f));

        std::cout << "tmp2 : " << std::endl;
        VecView(tmp1, PETSC_VIEWER_STDOUT_WORLD);
        std::cout << std::endl;

        std::cout << "A : " << std::endl;
        MatView(ls->A, PETSC_VIEWER_STDOUT_WORLD);
        std::cout << std::endl;

        KSPSolve(ksp,tmp2,tmp1); 

        std::cout << "tmp1 : " << std::endl;
        VecView(tmp1, PETSC_VIEWER_STDOUT_WORLD);
        std::cout << std::endl;

        PetscCall(MatMult(ls->A,tmp1,test));

        std::cout << "test : " << std::endl;
        VecView(test, PETSC_VIEWER_STDOUT_WORLD);
        std::cout << std::endl;

        PetscCall(VecAXPY(ls->x,1,tmp1)); // x1

        //std::cout << "tmp1 : " << std::endl;
        //VecView(tmp1, PETSC_VIEWER_STDOUT_WORLD);
        //std::cout << std::endl;

        //std::cout << "x : " << std::endl;
        //VecView(ls->x, PETSC_VIEWER_STDOUT_WORLD);
        //std::cout << std::endl;

        PetscCall(MatMult(B,ls->x,tmp3));
        PetscCall(VecAXPY(tmp3, -1, ls->g));
        PetscCall(VecScale(tmp3,-1.0));

        //std::cout << "tmp3 : " << std::endl;
        //VecView(tmp3, PETSC_VIEWER_STDOUT_WORLD);
        //std::cout << std::endl;

        PetscCall(VecAXPY(ls->y,tau,tmp3));

        //std::cout << "y : " << std::endl;
        //VecView(tmp3, PETSC_VIEWER_STDOUT_WORLD);
        //std::cout << std::endl;

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
