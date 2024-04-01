#include "solve.h"

PetscErrorCode CreateLinearSys(linearSys * ls, ReducedSys * reducedsys){

    // Compute target linear system
    // returns
    // A B F
    // B O G
    // C will be assigned separately
    int M, N;

    // Boundary term right hand side
    // velocity
    Vec g1;
    PetscCall(MatGetSize(reducedsys->M, &M, &N));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &g1));
    PetscCall(VecSetSizes(g1, PETSC_DECIDE, M));
    PetscCall(VecSetUp(g1));
  
    PetscCall(MatMult(reducedsys->Kg, reducedsys->g, g1));

    // pressure
    Vec g2;
    PetscCall(MatGetSize(reducedsys->B, &M, &N));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &g2));
    PetscCall(VecSetSizes(g2, PETSC_DECIDE, N));
    PetscCall(VecSetUp(g2));

    Mat BgT;
    PetscCall(MatCreateTranspose(reducedsys->Bg, &BgT));

    PetscCall(MatMult(BgT, reducedsys->g, g2));

    // Move boundary condition vectors to the right hand side of the 
    // rhs = g-dirichlet+neumann
    PetscCall(VecAYPX(g1,-1,reducedsys->source));

    //PetscCall(VecAXPY(g1,1.0,reducedsys->neum));

    PetscCall(VecScale(g2,-1));

    // Copy computed vectors to 
    PetscCall(MatConvert(reducedsys->B, MATSAME, MAT_INITIAL_MATRIX, &ls->B));
    PetscCall(MatConvert(reducedsys->M, MATSAME, MAT_INITIAL_MATRIX, &ls->A));

    PetscCall(VecDuplicate(g1,&ls->f));
    PetscCall(VecDuplicate(g2,&ls->g));

    VecCopy(g1, ls->f);
    VecCopy(g2, ls->g);

    PetscCall(VecCreate(PETSC_COMM_WORLD, &ls->x));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &ls->y));

    PetscCall(VecSetSizes(ls->x,PETSC_DECIDE,M));
    PetscCall(VecSetSizes(ls->y,PETSC_DECIDE,N));

    PetscCall(VecSetUp(ls->x));
    PetscCall(VecSetUp(ls->y));

    PetscCall(VecCopy(g1,ls->x));
    PetscCall(VecCopy(g2,ls->y));
   
    PetscCall(VecZeroEntries(ls->x));
    PetscCall(VecZeroEntries(ls->y));

    PetscCall(MatCreate(PETSC_COMM_WORLD, &ls->C));
    PetscCall(MatSetSizes(ls->C, PETSC_DECIDE, PETSC_DECIDE, M*N, M*N));
    PetscCall(MatSetUp(ls->C));

    PetscCall(MatAssemblyBegin(ls->C, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(ls->C, MAT_FINAL_ASSEMBLY));

    PetscCall(MatZeroEntries(ls->C));

    return PETSC_SUCCESS;
}

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

    while(r > tol && iter < MaxIter){

        PetscCall(MatMult(ls->B, ls->y, tmp1));
        PetscCall(MatMult(ls->A, ls->x, tmp2));

        PetscCall(VecAXPY(tmp2,-1.0,tmp1));

        PetscCall(VecAYPX(tmp2,-1.0,ls->f));

        KSPSolve(ksp,tmp2,tmp1); 

        PetscCall(VecAXPY(ls->x,1,tmp1)); // x1

        PetscCall(MatMult(B,ls->x,tmp3));
        PetscCall(VecAXPY(tmp3, -1, ls->g));
        PetscCall(VecScale(tmp3,-1.0));

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

PetscErrorCode InexactUzawa(linearSys * ls, double tol, int MaxIter, double tau1, double tau2){

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

    Vec tmp1, tmp2, tmp3, tmp4;

    PetscCall(VecDuplicate(ls->f, &tmp1));
    PetscCall(VecDuplicate(ls->f, &tmp2));
    PetscCall(VecDuplicate(ls->g, &tmp3));
    PetscCall(VecDuplicate(ls->g, &tmp4));

    PetscCall(VecZeroEntries(tmp1));
    PetscCall(VecZeroEntries(tmp2));
    PetscCall(VecZeroEntries(tmp3));
    PetscCall(VecZeroEntries(tmp4));

    int size;
    VecGetSize(tmp3,&size);

    double *arraytmp; 

    while(r > tol && iter < MaxIter){

        PetscCall(MatMult(ls->B, ls->y, tmp1));
        PetscCall(MatMult(ls->A, ls->x, tmp2));

        PetscCall(VecAXPY(tmp2,-1.0,tmp1));

        PetscCall(VecAYPX(tmp2,-1.0,ls->f));

        KSPSolve(ksp,tmp2,tmp1); 

        PetscCall(VecAXPY(ls->x,1,tmp1)); // x1

        PetscCall(MatMult(B,ls->x,tmp3));
        PetscCall(MatMult(ls->C,ls->y,tmp4));
        PetscCall(VecAXPY(tmp3,1.0,tmp4));
        PetscCall(VecAXPY(tmp3, -1, ls->g));
        PetscCall(VecScale(tmp3,-1.0));

/*
        // Take out tmp3 constant kernal =====
        double mean;
        VecMean(tmp3, &mean);
        PetscCall(VecGetArray(tmp3, &arraytmp));
        for (unsigned int k=0; k<size; k++){
            arraytmp[k] -= mean;
        }
        PetscCall(VecRestoreArray(tmp3, &arraytmp));

        // ===================================
*/
        // Uniform tau value
        //PetscCall(VecAXPY(ls->y,tau,tmp3));

        // Different tau values
        Vec tmp31, tmp32;
        PetscCall(VecNestGetSubVec(tmp3, 0, &tmp31));
        PetscCall(VecNestGetSubVec(tmp3, 1, &tmp32));

        PetscCall(VecScale(tmp31,tau1));
        PetscCall(VecScale(tmp32,tau2));

        //PetscCall(VecView(tmp3,PETSC_VIEWER_STDOUT_WORLD));

        PetscCall(VecAXPY(ls->y,1.0,tmp3));

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


PetscErrorCode CoupledSolver(linearSys * ls1, linearSys * ls2, linearSys * lsResult, Mat * K,
                             double tol, int MaxIter, double tau1, double tau2){

    // Solve coupled system with Uzawa algorithm
    // Create coupled system with two different linear system
    // The original system 
    // As  -Bs               
    // BsT  Cs       K      
    //         Ad   -Bd     
    //      K  BdT   Cd
    // Symmetrically permutate the original system 
    // As     | -Bs        xs     ubs
    // ____Ad_|_____-Bd  * xd  =  ubd 
    // BsT    |  Cs  K     ys     qs
    //     BdT|  K   Cd    yd     qd
    // Form a new saddle point system

    int M1, N1, M2, N2;

    PetscCall(VecGetSize(ls1->f, &M1)); 
    PetscCall(VecGetSize(ls1->g, &N1));
    PetscCall(VecGetSize(ls2->f, &M2)); 
    PetscCall(VecGetSize(ls2->g, &N2));

    // Create coupled A matrix 
    Mat arrayA[4], Z, ZT;

    // Create Two zero matrix
    PetscCall(MatCreate(PETSC_COMM_WORLD, &Z));
    PetscCall(MatSetSizes(Z,PETSC_DECIDE,PETSC_DECIDE,M1,M2));
    PetscCall(MatSetUp(Z));
    PetscCall(MatAssemblyBegin(Z,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Z,MAT_FINAL_ASSEMBLY));
    PetscCall(MatZeroEntries(Z));
    PetscCall(MatCreateTranspose(Z,&ZT));
   
    arrayA[0] = ls1->A;
    arrayA[1] = Z;
    arrayA[2] = ZT;
    arrayA[3] = ls2->A;

    PetscCall(MatCreateNest(PETSC_COMM_WORLD,2, NULL, 2, NULL, arrayA, &lsResult->A));

    // Create Coupled B matrix 
    Mat arrayB[4], Zb1, Zb2;

    PetscCall(MatCreate(PETSC_COMM_WORLD, &Zb1));
    PetscCall(MatSetSizes(Zb1,PETSC_DECIDE, PETSC_DECIDE, M1, N2));
    PetscCall(MatSetUp(Zb1));
    PetscCall(MatAssemblyBegin(Zb1,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Zb1,MAT_FINAL_ASSEMBLY));
    PetscCall(MatZeroEntries(Zb1));

    PetscCall(MatCreate(PETSC_COMM_WORLD, &Zb2));
    PetscCall(MatSetSizes(Zb2,PETSC_DECIDE, PETSC_DECIDE, M2, N1));
    PetscCall(MatSetUp(Zb2));
    PetscCall(MatAssemblyBegin(Zb2,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Zb2,MAT_FINAL_ASSEMBLY));
    PetscCall(MatZeroEntries(Zb2));

    arrayB[0] = ls1->B;
    arrayB[1] = Zb1;
    arrayB[2] = Zb2;
    arrayB[3] = ls2->B;

    PetscCall(MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,arrayB, &lsResult->B));

    // Create Coupled C matrix
    Mat arrayC[4];

    arrayC[0] = ls1->C;
    arrayC[1] = *K;
    arrayC[2] = *K;
    arrayC[3] = ls2->C;

    PetscCall(MatCreateNest(PETSC_COMM_WORLD,2, NULL, 2, NULL, arrayC, &lsResult->C));

    Vec arrayx[2], arrayy[2], arrayf[2], arrayg[2];

    arrayx[0] = ls1->x;
    arrayx[1] = ls2->x;

    arrayy[0] = ls1->y;
    arrayy[1] = ls2->y;

    arrayf[0] = ls1->f;
    arrayf[1] = ls2->f;

    arrayg[0] = ls1->g;
    arrayg[1] = ls2->g;

    PetscCall(VecCreateNest(PETSC_COMM_WORLD,2,NULL,arrayx,&lsResult->x));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD,2,NULL,arrayy,&lsResult->y));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD,2,NULL,arrayf,&lsResult->f));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD,2,NULL,arrayg,&lsResult->g));

    return InexactUzawa(lsResult, tol, MaxIter, tau1, tau2);
}
