#include "psolve.h"

int CreateLinearSys(ReducedSys * redsys, const int& nelem){

    // Right hand side F
    PetscCall(VecDuplicate(redsys->source, &redsys->F));
    PetscCall(MatMult(redsys->Kg, redsys->g, redsys->F));

    PetscCall(VecAYPX(redsys->F, -1, redsys->source));

    // Right hand side G
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nelem, &redsys->G));

    Mat BgT;
    PetscCall(MatCreateTranspose(redsys->Bg, &BgT));
    PetscCall(MatMult(BgT, redsys->g, redsys->G));
    PetscCall(VecScale(redsys->G, -1));

    return 0;
}

int CreateCoupledSystem(ReducedSys * redsys1, ReducedSys * redsys2,
                        ReducedSys * redsysRes, Mat * K){


    // Create a coupled saddle point system
    // with two small saddle point system
    int M1, N1, M2, N2;
    PetscCall(VecGetSize(redsys1->F, &M1)); 
    PetscCall(VecGetSize(redsys1->G, &N1));
    PetscCall(VecGetSize(redsys2->F, &M2)); 
    PetscCall(VecGetSize(redsys2->G, &N2));

    // Create a Coupled A matrix
    Mat arrayA[4], Z, ZT;

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                           M1, M2, 0, NULL, 0, NULL, &Z));
    PetscCall(MatSetUp(Z));

    PetscCall(MatAssemblyBegin(Z,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Z,MAT_FINAL_ASSEMBLY));
 
    PetscCall(MatCreateTranspose(Z,&ZT));
 
    arrayA[0] = redsys1->M;
    arrayA[1] = Z;
    arrayA[2] = ZT;
    arrayA[3] = redsys2->M;

    PetscCall(MatCreateNest(PETSC_COMM_WORLD,2, NULL, 2, NULL, arrayA, 
                            &redsysRes->M));

    // Create Coupled B matrix 
    Mat arrayB[4], Zb1, Zb2;
    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                           M1, N2, 0, NULL, 0, NULL, &Zb1));
    PetscCall(MatSetUp(Zb1));

    PetscCall(MatAssemblyBegin(Zb1,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Zb1,MAT_FINAL_ASSEMBLY));
 
    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                           M2, N1, 0, NULL, 0, NULL, &Zb2));
    PetscCall(MatSetUp(Zb2));

    PetscCall(MatAssemblyBegin(Zb2,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Zb2,MAT_FINAL_ASSEMBLY));
 
    arrayB[0] = redsys1->B;
    arrayB[1] = Zb1;
    arrayB[2] = Zb2;
    arrayB[3] = redsys2->B;

    PetscCall(MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,arrayB, &redsysRes->B));

    // Create Coupled C matrix
    Mat arrayC[4];

    arrayC[0] = redsys1->C;
    arrayC[1] = *K;
    arrayC[2] = *K;
    arrayC[3] = redsys2->C;

    PetscCall(MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,arrayC, &redsysRes->C));

    // Create right hand side and solution nested vectors
    Vec arrayx[2], arrayy[2], arrayf[2], arrayg[2];

    arrayf[0] = redsys1->F;
    arrayf[1] = redsys2->F;

    arrayg[0] = redsys1->G;
    arrayg[1] = redsys2->G;

    PetscCall(VecCreateNest(PETSC_COMM_WORLD,2,NULL,arrayf,&redsysRes->F));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD,2,NULL,arrayg,&redsysRes->G));

    PetscCall(VecDuplicate(redsysRes->F, &redsysRes->x));
    PetscCall(VecDuplicate(redsysRes->G, &redsysRes->y));

    return 0;
}

int CoupledUzawa(ReducedSys * redsys, double tol, int MaxIter){
    // final form of the coupled uzawa iteration
    // no need of tau1 and tau2 parameter

    MatScale(redsys->B, -1);
    MatScale(redsys->C, -1);
    VecScale(redsys->G, -1);

    // ===========================================================
    KSP kspCG;
    PC  pcCG; 
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspCG));
    PetscCall(KSPSetOperators(kspCG, redsys->M, redsys->M));
    PetscCall(KSPSetType(kspCG, KSPCG));
    PetscCall(KSPCGSetType(kspCG, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspCG, PETSC_FALSE));
 
    // Create B transpose
    Mat BT;
    PetscCall(MatCreateTranspose(redsys->B, &BT));

    // Define KSP for schur complement for Darcy part
    KSP kspMINRESd, kspSchurd, kspMINRESs, kspSchurs;
    Mat Sd, Ad, Bd, BdT, Cd;
    Mat Ss, As, Bs, BsT, Cs;

    PetscCall(MatNestGetSubMat(redsys->M, 1, 1, &Ad));
    PetscCall(MatNestGetSubMat(redsys->B, 1, 1, &Bd));
    PetscCall(MatNestGetSubMat(redsys->C, 1, 1, &Cd));
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
    PetscCall(KSPSetTolerances(kspMINRESd, 1e-25, 10e-20, 10, 2000));

    // ===================================================================
    PetscCall(MatNestGetSubMat(redsys->M, 0, 0, &As));
    PetscCall(MatNestGetSubMat(redsys->B, 0, 0, &Bs));
    PetscCall(MatNestGetSubMat(redsys->C, 0, 0, &Cs));
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
    PetscCall(KSPSetTolerances(kspMINRESs, 1e-25, 10e-20, 10, 2000));

    double r = 1.0;
    int    iter = 0;

    Vec tmp1, tmp2, tmp3, tmp4;

    PetscCall(VecDuplicate(redsys->F, &tmp1));
    PetscCall(VecDuplicate(redsys->F, &tmp2));
    PetscCall(VecDuplicate(redsys->G, &tmp3));
    PetscCall(VecDuplicate(redsys->G, &tmp4));

    PetscCall(VecZeroEntries(tmp1));
    PetscCall(VecZeroEntries(tmp2));
    PetscCall(VecZeroEntries(tmp3));
    PetscCall(VecZeroEntries(tmp4));

    PetscCall(VecZeroEntries(redsys->x));
    PetscCall(VecZeroEntries(redsys->y));

    Vec tmp31, tmp32;

    // Iteration starts here

    while(r>tol && iter < MaxIter){
 
        PetscCall(MatMult(redsys->B, redsys->y, tmp1));

        PetscCall(MatMult(redsys->M, redsys->x, tmp2));

        // tmp2 = ls-f - (Ax + B'y)
        PetscCall(VecAXPBYPCZ(tmp2, 1.0, -1.0, -1.0, redsys->F, tmp1));

        // tmp1 = A^-1 tmp2
        PetscCall(KSPSolve(kspCG,tmp2,tmp1));

        PetscCall(VecAXPY(redsys->x,1,tmp1)); 

        PetscCall(MatMult(BT,redsys->x,tmp3));
        PetscCall(MatMult(redsys->C,redsys->y,tmp4));

        // tmp3 = Bx + Cy - G
        PetscCall(VecAXPBYPCZ(tmp3, -1.0, 1.0, 1.0, redsys->G, tmp4));

        // Use MINRES to calculate Darcy part 
        PetscCall(VecNestGetSubVec(tmp3, 0, &tmp31));
        PetscCall(VecNestGetSubVec(tmp3, 1, &tmp32));
                  
        KSPSolve(kspMINRESs, tmp31, tmp31);
        KSPSolve(kspMINRESd, tmp32, tmp32);

        PetscCall(VecAXPY(redsys->y,-1.0,tmp3));

        // Check norm of increment
        PetscReal val1, val2;
        PetscCall(VecNorm(tmp1,NORM_2,&val1));
        PetscCall(VecNorm(tmp3,NORM_2,&val2));
        r = val1 + val2; 

        iter++;

    }

    Vec tmpDarcy;
    PetscCall(VecNestGetSubVec(redsys->x, 1, &tmpDarcy));
    VecScale(tmpDarcy, -1);  

    if (iter < MaxIter){
        printf("Uzawa converged successfully! r = %.3e, Used %d iterations. \n", r, iter);
        return 0;
    } else {
        printf("Uzawa failed to converge! r = %.3e \n", r);
        return -1;
    }

}

int Uzawa(ReducedSys * redsys, double tol, int MaxIter){

    MatScale(redsys->B, -1);
    MatScale(redsys->C, -1);
    VecScale(redsys->G, -1);
//MatView(redsys->C, PETSC_VIEWER_STDOUT_WORLD);
//    MatZeroEntries(redsys->C);

    // ===========================================================
    KSP kspCG;
    PC  pcCG; 
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspCG));
    PetscCall(KSPSetOperators(kspCG, redsys->M, redsys->M));
    PetscCall(KSPSetType(kspCG, KSPCG));
    PetscCall(KSPCGSetType(kspCG, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspCG, PETSC_FALSE));

    // Define KSP for schur complement for Darcy part
    KSP kspMINRES, kspSchur;
    Mat S, BT;

    PetscCall(MatCreateTranspose(redsys->B, &BT));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspMINRES));
    PetscCall(MatCreateSchurComplement(redsys->M, redsys->M, redsys->B, BT, redsys->C, &S));
 
    PetscCall(MatSchurComplementGetKSP(S, &kspSchur));
    PetscCall(KSPSetType(kspSchur, KSPCG));
    PetscCall(KSPCGSetType(kspSchur, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspSchur, PETSC_FALSE));

    PetscCall(KSPSetOperators(kspMINRES, S, S));
    PetscCall(KSPSetType(kspMINRES, KSPMINRES)); 
    PetscCall(KSPSetInitialGuessNonzero(kspMINRES, PETSC_FALSE));
    PetscCall(KSPSetTolerances(kspMINRES, 1e-16, 10e-20, 10, 2000));

    double r = 1.0;
    int    iter = 0;

    Vec tmp1, tmp2, tmp3, tmp4;

    PetscCall(VecDuplicate(redsys->F, &tmp1));
    PetscCall(VecDuplicate(redsys->F, &tmp2));
    PetscCall(VecDuplicate(redsys->G, &tmp3));
    PetscCall(VecDuplicate(redsys->G, &tmp4));

    PetscCall(VecZeroEntries(tmp1));
    PetscCall(VecZeroEntries(tmp2));
    PetscCall(VecZeroEntries(tmp3));
    PetscCall(VecZeroEntries(tmp4));

    PetscCall(VecDuplicate(redsys->F, &redsys->x));
    PetscCall(VecDuplicate(redsys->G, &redsys->y));

    PetscCall(VecZeroEntries(redsys->x));
    PetscCall(VecZeroEntries(redsys->y));

    while (r> tol && iter < MaxIter){

        PetscCall(MatMult(redsys->B, redsys->y, tmp1));

        PetscCall(MatMult(redsys->M, redsys->x, tmp2));

        // tmp2 = ls-f - (Ax + B'y)
        PetscCall(VecAXPBYPCZ(tmp2, 1.0, -1.0, -1.0, redsys->F, tmp1));

        // tmp1 = A^-1 tmp2
        PetscCall(KSPSolve(kspCG,tmp2,tmp1));

        PetscCall(VecAXPY(redsys->x,1,tmp1)); 

        PetscCall(MatMult(BT,redsys->x,tmp3));
        PetscCall(MatMult(redsys->C,redsys->y,tmp4));

        // tmp3 = Bx + Cy - G
        PetscCall(VecAXPBYPCZ(tmp3, -1.0, 1.0, 1.0, redsys->G, tmp4));

        // Use MINRES to calculate Darcy part 
        KSPSolve(kspMINRES, tmp3, tmp3);

        PetscCall(VecAXPY(redsys->y,-1.0,tmp3));

        // Check norm of increment
        PetscReal val1, val2;
        PetscCall(VecNorm(tmp1,NORM_2,&val1));
        PetscCall(VecNorm(tmp3,NORM_2,&val2));
        r = val1 + val2; 

        iter++;

    }

    if (iter < MaxIter){
        printf("Uzawa converged successfully! r = %.3e, Used %d iterations. \n", r, iter);
        return 0;
    } else {
        printf("Uzawa failed to converge! r = %.3e \n", r);
        return -1;
    }
}
