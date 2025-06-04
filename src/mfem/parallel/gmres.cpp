#include "psolve.h"

int SchurSolver(ReducedSys * redsys){

    // Create a big nested matrix and feed it directly into GMRES
    // with two small saddle point system

    MatScale(redsys->B, -1);
    MatScale(redsys->C, -1);
    VecScale(redsys->G, -1);

    // Create B transpose
    Mat BT;
    PetscCall(MatCreateTranspose(redsys->B, &BT));

    KSP kspA;
    PC  pcA;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspA));
    PetscCall(KSPSetOperators(kspA, redsys->M, redsys->M));
    PetscCall(KSPSetType(kspA, KSPCG));
    PetscCall(KSPCGSetType(kspA, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspA, PETSC_FALSE));
    PetscCall(KSPSetTolerances(kspA, 1e-25, 10e-20, 10, 2000));

    // Schur Complement S = (C - BTA^(-1)B)
    Mat S;
    KSP kspSchur;
    PetscCall(MatCreateSchurComplement(redsys->M, redsys->M, redsys->B, BT, redsys->C, &S));
    PetscCall(MatSchurComplementGetKSP(S, &kspSchur));
    PetscCall(KSPSetType(kspSchur, KSPCG));
    PetscCall(KSPCGSetType(kspSchur, KSP_CG_SYMMETRIC));
    PetscCall(KSPSetInitialGuessNonzero(kspSchur, PETSC_FALSE));
    PetscCall(KSPSetTolerances(kspSchur, 1e-25, 1e-20, 10, 2000));

    
    // Compute g-BTA^(-1)f
    Vec tmp1, tmp2, tmp3;
    PetscCall(VecDuplicate(redsys->F, &tmp1));
    PetscCall(VecDuplicate(redsys->F, &tmp2));
    PetscCall(VecDuplicate(redsys->G, &tmp3));
    PetscCall(VecCopy(redsys->G, tmp3));
 
    PetscCall(KSPSolve(kspA, redsys->F, tmp1));
    PetscCall(MatMult(BT, tmp1, tmp2));
    PetscCall(VecAXPY(tmp3, -1, tmp2));

    // Solver for y = S^(-1)
    KSP ksp2;
    PetscCall(KSPSetOperators(ksp2, S, S));
    PetscCall(KSPSetType(ksp2, KSPMINRES));

    


    return 0;
}
