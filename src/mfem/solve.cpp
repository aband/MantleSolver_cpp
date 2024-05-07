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

    // Divergence free option
    PetscCall(VecScale(g2,-1));

    // Divergence not free
    //PetscCall(VecAYPX(g2,-1,reducedsys->source))

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

PetscErrorCode SimpleUzawa(linearSys * ls, 
                           double tol, int MaxIter, 
                           double tau, int precondType){

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

        switch(precondType) {
            case 0:
                // Choice 1
                // Uniform tau value
                PetscCall(VecScale(tmp3, tau));
                break;

            case 1:
                KSP kspQ;
                Mat BT, Sd; 

                PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspQ));
                PetscCall(MatCreateTranspose(ls->B, &BT));

                PetscCall(MatCreateSchurComplement(ls->A, ls->A, ls->B, BT, NULL , &Sd));

                KSP kspschur;
                PC  pc;
                PetscCall(MatSchurComplementGetKSP(Sd, &kspschur));
                KSPSetType(kspschur,KSPCG);
                KSPCGSetType(kspschur,KSP_CG_SYMMETRIC);
                KSPSetInitialGuessNonzero(kspschur, PETSC_FALSE);
                KSPSetTolerances(kspschur, 10e-7, 10e-12, 5, 2);

                //MatView(Sd, PETSC_VIEWER_STDOUT_WORLD); 

                PetscCall(KSPSetOperators(kspQ, Sd, Sd));
                PetscCall(KSPSetType(kspQ, KSPGMRES));
                PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE));

                VecView(tmp3, PETSC_VIEWER_STDOUT_WORLD);

                KSPSolve(kspQ, tmp3, tmp3);

                VecView(tmp3, PETSC_VIEWER_STDOUT_WORLD);

                break;

        }

//        PetscCall(VecAXPY(ls->y,tau,tmp3));
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

PetscErrorCode InexactUzawa(linearSys * ls, double tol, int MaxIter, double tau1, double tau2, int flag){

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
    //KSPSetTolerances(ksp, 10e-7, 10e-7, 10e-7, 10000);
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

        VecView(tmp3, PETSC_VIEWER_STDOUT_WORLD);
        VecView(tmp4, PETSC_VIEWER_STDOUT_WORLD);

        PetscCall(VecAXPY(tmp3,1.0,tmp4));
        PetscCall(VecAXPY(tmp3, -1, ls->g));
        PetscCall(VecScale(tmp3,-1.0));

/* ====================================================================================
        // Take out tmp3 constant kernal =====
        // Should not be used if uzawa iteration is implemented correctly
        double mean;
        VecMean(tmp3, &mean);
        PetscCall(VecGetArray(tmp3, &arraytmp));
        for (unsigned int k=0; k<size; k++){
            arraytmp[k] -= mean;
        }
        PetscCall(VecRestoreArray(tmp3, &arraytmp));
======================================================================================= */

        Vec tmp31, tmp32;

        KSP kspMINRES, kspSchur;

        switch(flag){
            case 0:
                // Choice 1
                // Uniform tau value
                PetscCall(VecScale(tmp3, tau1));
                break;

            case 1:
                // Choice 2
                // Different tau values
 
                PetscCall(VecNestGetSubVec(tmp3, 0, &tmp31));
                PetscCall(VecNestGetSubVec(tmp3, 1, &tmp32));

                PetscCall(VecScale(tmp31,tau1));
                PetscCall(VecScale(tmp32,tau2));

                //PetscCall(VecView(tmp3,PETSC_VIEWER_STDOUT_WORLD));
                break;

            case 2:
                // Choice 3
                KSP kspQ;
                Mat AQ, As, Ad, Bd, BdT, tmp, Sd, Cd;
                int M, N;
                Vec D, xd;

                PetscCall(VecNestGetSubVec(ls->x,1,&xd));
                PetscCall(VecDuplicate(xd, &D));

                // Extracting diagonal of A as approximation of A matrix
                PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspQ));
                // Extract Ad and Bd from sub matrix
                PetscCall(MatNestGetSubMat(ls->A, 1, 1, &Ad));
                PetscCall(MatNestGetSubMat(ls->B, 1, 1, &Bd));
                PetscCall(MatNestGetSubMat(ls->C, 1, 1, &Cd));
                PetscCall(MatCreateTranspose(Bd, &BdT));

                PetscCall(MatGetSize(Ad, &M, &N));

                // Extract diagonal from Darcy mass matrix
                // Get the inverse value from it
                PetscCall(MatGetDiagonal(Ad, D));
                PetscCall(VecReciprocal(D));

                PetscCall(MatDuplicate(Ad, MAT_DO_NOT_COPY_VALUES, &AQ));

                PetscCall(MatDiagonalSet(AQ, D, INSERT_VALUES));

                // Compute schur complement explicitly
                // Sd = BT A^-1 B + C
                PetscCall(MatMatMult(AQ, Bd, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmp));

                PetscCall(MatMatMult(BdT, tmp, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Sd));
                PetscCall(MatAXPY(Sd, 1.0, Cd, DIFFERENT_NONZERO_PATTERN));

                PetscCall(VecNestGetSubVec(tmp3, 0, &tmp31));
                PetscCall(VecNestGetSubVec(tmp3, 1, &tmp32));
                //PetscCall(KSPSetOperators(kspQ, ))
                PetscCall(VecScale(tmp31, tau1));

                //MatView(Sd, PETSC_VIEWER_STDOUT_WORLD);

                PetscCall(KSPSetOperators(kspQ, Sd, Sd));
                PetscCall(KSPSetType(kspQ, KSPMINRES));
                PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE));

                KSPSolve(kspQ, tmp32, tmp32);

                break;

            default:
                cout << "No valid flad defined in Uzawa Solver." << endl;
                break;
        }


        PetscCall(VecAXPY(ls->y,1.0,tmp3));

        // Check norm of increment
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

PetscErrorCode InexactUzawa(linearSys * ls, double tol, int MaxIter, double tau){

    // Create Schur complement for Darcy part

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
    //PetscCall(KSPSetTolerances(ksp,10e-7, 10e-12, 10e-12, 100));

    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCBJACOBI));

    Mat B;
    PetscCall(MatCreateTranspose(ls->B, &B));

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
    VecGetSize(tmp3, &size);

    // Uzawa loop
    while(r>tol && iter < MaxIter){
        PetscCall(MatMult(ls->B, ls->y, tmp1));
        PetscCall(MatMult(ls->A, ls->x, tmp2));

        PetscCall(VecAXPY(tmp2,-1.0,tmp1));

        PetscCall(VecAYPX(tmp2,-1.0,ls->f));

        KSPSolve(ksp,tmp2,tmp1); 

        VecView(tmp1, PETSC_VIEWER_STDOUT_WORLD);

        PetscCall(VecAXPY(ls->x,1,tmp1)); // x1

        PetscCall(MatMult(B,ls->x,tmp3));
        PetscCall(MatMult(ls->C,ls->y,tmp4));

        //VecView(tmp3, PETSC_VIEWER_STDOUT_WORLD);
        //VecView(tmp4, PETSC_VIEWER_STDOUT_WORLD);
        //VecView(ls->g, PETSC_VIEWER_STDOUT_WORLD);

        PetscCall(VecAXPY(tmp3,1.0,tmp4));
        PetscCall(VecAXPY(tmp3, -1, ls->g));
        PetscCall(VecScale(tmp3,-1.0));

        PetscCall(VecAXPY(ls->y,tau,tmp3));

        PetscReal val1, val2;
        PetscCall(VecNorm(tmp1,NORM_2,&val1));
        PetscCall(VecNorm(tmp3,NORM_2,&val2));
        r = val1 + val2; 

        //cout << val1 << "  " << val2 << endl << endl;

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
                             double tol, int MaxIter, double tau1, double tau2, int flag){

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

    //return InexactUzawa(lsResult, tol, MaxIter, tau1, tau2, flag);
    //return ExactUzawa(lsResult, tol, MaxIter);
    return CoupledExactUzawa(lsResult, tau1, tau2, tol, MaxIter, flag);
}

PetscErrorCode CoupledSolver(linearSys * ls1, linearSys * ls2, linearSys * lsResult, Mat * K, 
                             double tol, int MaxIter, double tau){

    // Convert Darcy part to corresponding schur complement
    // The new resulting system
    // As  -Bs        xs    fs
    // BsT  Cs  K   * ys  = gs
    //      K   Sd    yd    gd - BdT Ad^-1 fd
    // Form a smaller saddle point system
    // Ad is guaranteed positive definite
    
    int M1, N1, M2, N2;
    PetscCall(VecGetSize(ls1->f,&M1));
    PetscCall(VecGetSize(ls1->g,&N1));
    PetscCall(VecGetSize(ls2->f,&M2));
    PetscCall(VecGetSize(ls2->g,&N2));

    // Create A matrix for linear system
    PetscCall(MatConvert(ls1->A, MATSAME, MAT_INITIAL_MATRIX, &lsResult->A));

    // Create B matrix for linear system 
    Mat arrayB[2], Zb; 
    PetscCall(MatCreate(PETSC_COMM_WORLD, &Zb));
    PetscCall(MatSetSizes(Zb,PETSC_DECIDE, PETSC_DECIDE, M1, N2));
    PetscCall(MatSetUp(Zb));
    PetscCall(MatAssemblyBegin(Zb,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Zb,MAT_FINAL_ASSEMBLY));
    PetscCall(MatZeroEntries(Zb));

    arrayB[0] = ls1->B;
    arrayB[1] = Zb;

    PetscCall(MatCreateNest(PETSC_COMM_WORLD, 1, NULL, 2, NULL, arrayB, &lsResult->B));

    // Create Coupled C matrix
    Mat arrayC[4];

    // Create schur complement for Darcy system
    Mat Sd,BT,mC;

    PetscCall(MatCreateTranspose(ls2->B, &BT));

    PetscCall(MatDuplicate(ls2->C, MAT_COPY_VALUES, &mC));
    PetscCall(MatScale(mC, -1));

    PetscCall(MatCreateSchurComplement(ls2->A, ls2->A, ls2->B, BT, mC, &Sd));

    KSP kspschur;
    PC  pc;
    PetscCall(MatSchurComplementGetKSP(Sd, &kspschur));
    KSPSetType(kspschur,KSPCG);
    KSPCGSetType(kspschur,KSP_CG_SYMMETRIC);
    KSPSetInitialGuessNonzero(kspschur, PETSC_FALSE);
    KSPSetTolerances(kspschur, 10e-7, 10e-12, 5, 1000);

    PetscCall(KSPGetPC(kspschur, &pc));
    PetscCall(PCSetType(pc, PCBJACOBI));

    MatView(Sd, PETSC_VIEWER_STDOUT_WORLD);

    arrayC[0] = ls1->C;
    arrayC[1] = *K;
    arrayC[2] = *K;
    arrayC[3] = Sd;

    PetscCall(MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, arrayC, &lsResult->C));

    // Create right hand side vectors
    // new right hand side gd - BdT Ad-1 fd

    Vec arrayy[2];

    PetscCall(VecDuplicate(ls1->f, &lsResult->f));
    PetscCall(VecCopy(ls1->f, lsResult->f));

    arrayy[0] = ls1->g;
    arrayy[1] = ls2->g; // temperory
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, arrayy, &lsResult->g));

    // Create solution vectors
    PetscCall(VecDuplicate(ls1->x, &lsResult->x));
    PetscCall(VecZeroEntries(lsResult->x));

    PetscCall(VecDuplicate(lsResult->g, &lsResult->y));
    PetscCall(VecZeroEntries(lsResult->y));

    //return InexactUzawa(lsResult, tol, MaxIter, tau);
    return ExactUzawa(lsResult, tol, MaxIter);
}

PetscErrorCode RetrieveDarcy(linearSys * ls, Vec * darcy){

    // Retrieve darcy velocity from schur complement

    KSP ksp;
    PC  pc;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, ls->A, ls->A));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPCGSetType(ksp, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE)); // zero initial guess
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCBJACOBI));

    Vec darcyy;

    PetscCall(VecNestGetSubVec(ls->y,1,&darcyy));


    return PETSC_SUCCESS;
}

// ================================================================================

PetscErrorCode ExactUzawa(linearSys * ls, double tol, int MaxIter){
    // Variant of traditional uzawa solver
    // Solve singular schur complement system with MINRES
   
    MatScale(ls->B, -1);
    MatScale(ls->C, -1);
    VecScale(ls->g, -1);

    KSP kspCG;
    PC  pcCG; 
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspCG));
    PetscCall(KSPSetOperators(kspCG, ls->A, ls->A));
    PetscCall(KSPSetType(kspCG, KSPCG));
    PetscCall(KSPCGSetType(kspCG, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspCG, PETSC_FALSE));
    // Set up preconditioner
    //PetscCall(KSPSetTolerances(kspCG, 10e-7, 10e-7, 1000, 1000));
    //PetscCall(KSPGetPC(kspCG, &pcCG));
    //PetscCall(PCSetType(pcCG, PCBJACOBI));

    // Create B transpose
    Mat BT;
    PetscCall(MatCreateTranspose(ls->B, &BT));

    KSP  kspMINRES, kspSchur;
    PC   pcMINRES;
    Mat  S;

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspMINRES));
    PetscCall(MatCreateSchurComplement(ls->A, ls->A, ls->B, BT, ls->C, &S));
    
    PetscCall(MatSchurComplementGetKSP(S, &kspSchur));
    PetscCall(KSPSetType(kspSchur, KSPCG));
    PetscCall(KSPCGSetType(kspSchur, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspSchur, PETSC_FALSE));

    PetscCall(KSPSetOperators(kspMINRES, S, S));
    PetscCall(KSPSetType(kspMINRES, KSPMINRES)); 
    PetscCall(KSPSetInitialGuessNonzero(kspMINRES, PETSC_FALSE));
    PetscCall(KSPSetTolerances(kspMINRES, 1e-10, 1e-15, 1, 100));
 
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

    //MatView(ls->A, PETSC_VIEWER_STDOUT_WORLD);
    //MatView(ls->B, PETSC_VIEWER_STDOUT_WORLD);

    while(r > tol && iter < MaxIter){
        PetscCall(MatMult(ls->B, ls->y, tmp1));
        PetscCall(MatMult(ls->A, ls->x, tmp2));

        // tmp2 = ls-f - (Ax + B'y)
        PetscCall(VecAXPBYPCZ(tmp2, 1.0, -1.0, -1.0, ls->f, tmp1));

        // tmp1 = A^-1 tmp2
        PetscCall(KSPSolve(kspCG,tmp2,tmp1));
         
        PetscCall(VecAXPY(ls->x,1,tmp1)); 

        PetscCall(MatMult(BT,ls->x,tmp3));
        PetscCall(MatMult(ls->C,ls->y,tmp4));

        // tmp3 = Bx + Cy - G
        PetscCall(VecAXPBYPCZ(tmp3, -1.0, 1.0, 1.0, ls->g, tmp4));

        PetscCall(KSPSolve(kspMINRES,tmp3, tmp3));

        //KSPView(kspMINRES,PETSC_VIEWER_STDOUT_WORLD);

        PetscCall(VecAXPY(ls->y, -1.0, tmp3));

        // Check norm of increment
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

// Final form
PetscErrorCode CoupledExactUzawa(linearSys * ls, double tau1, double tau2,
                                 double tol, int MaxIter, int pType){

    // Exact Uzawa iteration for coupled system
    MatScale(ls->B, -1);
    MatScale(ls->C, -1);
    VecScale(ls->g, -1);

    // ===========================================================
    KSP kspCG;
    PC  pcCG; 
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspCG));
    PetscCall(KSPSetOperators(kspCG, ls->A, ls->A));
    PetscCall(KSPSetType(kspCG, KSPCG));
    PetscCall(KSPCGSetType(kspCG, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspCG, PETSC_FALSE));
    // Create B transpose
    Mat BT;
    PetscCall(MatCreateTranspose(ls->B, &BT));

    // Define KSP for schur complement for Darcy part
    KSP kspMINRESd, kspSchurd, kspMINRESs, kspSchurs;
    Mat Sd, Ad, Bd, BdT, Cd;
    Mat Ss, As, Bs, BsT, Cs;

    PetscCall(MatNestGetSubMat(ls->A, 1, 1, &Ad));
    PetscCall(MatNestGetSubMat(ls->B, 1, 1, &Bd));
    PetscCall(MatNestGetSubMat(ls->C, 1, 1, &Cd));
    PetscCall(MatCreateTranspose(Bd, &BdT));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspMINRESd));
    PetscCall(MatCreateSchurComplement(Ad, Ad, Bd, BdT, Cd, &Sd));
    
    PetscCall(MatSchurComplementGetKSP(Sd, &kspSchurd));
    PetscCall(KSPSetType(kspSchurd, KSPCG));
    PetscCall(KSPCGSetType(kspSchurd, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspSchurd, PETSC_FALSE));

    PetscCall(KSPSetOperators(kspMINRESd, Sd, Sd));
    PetscCall(KSPSetType(kspMINRESd, KSPMINRES)); 
    PetscCall(KSPSetInitialGuessNonzero(kspMINRESd, PETSC_FALSE));
    PetscCall(KSPSetTolerances(kspMINRESd, 1e-7, 10e-16, 10, 500));

    // ===================================================================
    PetscCall(MatNestGetSubMat(ls->A, 0, 0, &As));
    PetscCall(MatNestGetSubMat(ls->B, 0, 0, &Bs));
    PetscCall(MatNestGetSubMat(ls->C, 0, 0, &Cs));
    PetscCall(MatCreateTranspose(Bs, &BsT));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspMINRESs));
    PetscCall(MatCreateSchurComplement(As, As, Bs, BsT, Cs, &Ss));
    
    PetscCall(MatSchurComplementGetKSP(Ss, &kspSchurs));
    PetscCall(KSPSetType(kspSchurs, KSPCG));
    PetscCall(KSPCGSetType(kspSchurs, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspSchurs, PETSC_FALSE));

    PetscCall(KSPSetOperators(kspMINRESs, Ss, Ss));
    PetscCall(KSPSetType(kspMINRESs, KSPMINRES)); 
    PetscCall(KSPSetInitialGuessNonzero(kspMINRESs, PETSC_FALSE));
    PetscCall(KSPSetTolerances(kspMINRESs, 1e-7, 10e-16, 10, 500));

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

    Vec tmp31, tmp32;

    while(r>tol && iter < MaxIter){
        PetscCall(MatMult(ls->B, ls->y, tmp1));

        PetscCall(MatMult(ls->A, ls->x, tmp2));

        // tmp2 = ls-f - (Ax + B'y)
        PetscCall(VecAXPBYPCZ(tmp2, 1.0, -1.0, -1.0, ls->f, tmp1));

        // tmp1 = A^-1 tmp2
        PetscCall(KSPSolve(kspCG,tmp2,tmp1));

        PetscCall(VecAXPY(ls->x,1,tmp1)); 

        PetscCall(MatMult(BT,ls->x,tmp3));
        PetscCall(MatMult(ls->C,ls->y,tmp4));

        // tmp3 = Bx + Cy - G
        PetscCall(VecAXPBYPCZ(tmp3, -1.0, 1.0, 1.0, ls->g, tmp4));

        // Three different types of calculating increment vector
        switch(pType) {
            case 0:
                // A uniform diagonal style preconditioner
                PetscCall(VecScale(tmp3, tau1));

                break;

            case 1:
                // A diagonal style preoconditioner
                // but different values of tau for Stokes and Darcy problem

                PetscCall(VecNestGetSubVec(tmp3, 0, &tmp31));
                PetscCall(VecNestGetSubVec(tmp3, 1, &tmp32));

                PetscCall(VecScale(tmp31,tau1));
                PetscCall(VecScale(tmp32,tau2));

                break;

            case 2:

                // Use MINRES to calculate Darcy part 
                PetscCall(VecNestGetSubVec(tmp3, 0, &tmp31));
                PetscCall(VecNestGetSubVec(tmp3, 1, &tmp32));
                  
                KSPSolve(kspMINRESs, tmp31, tmp31);
                KSPSolve(kspMINRESd, tmp32, tmp32);

                break;

            default:

                PetscCall(PetscPrintf(PETSC_COMM_WORLD, "No valid flad defined in Uzawa Solver."));
 
                break;
        }

        PetscCall(VecAXPY(ls->y,-1.0,tmp3));

        // Check norm of increment
        PetscReal val1, val2;
        PetscCall(VecNorm(tmp1,NORM_2,&val1));
        PetscCall(VecNorm(tmp3,NORM_2,&val2));
        r = val1 + val2; 

        iter++;

    }

    Vec tmpDarcy;
    PetscCall(VecNestGetSubVec(ls->x, 1, &tmpDarcy));
    VecScale(tmpDarcy, -1);  

    if (iter < MaxIter){
        printf("Uzawa converged successfully! r = %.3e, Used %d iterations. \n", r, iter);
        return PETSC_SUCCESS;
    } else {
        printf("Uzawa failed to converge! r = %.3e \n", r);
        return PETSC_ERR_CONV_FAILED;
    }

}
