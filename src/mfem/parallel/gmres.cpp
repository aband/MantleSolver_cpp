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
    Vec tmp1, tmp2;
    PetscCall(VecDuplicate(redsys->F, &tmp1));
    PetscCall(VecDuplicate(redsys->G, &tmp2));
    PetscCall(VecDuplicate(redsys->G, &redsys->y));
    PetscCall(VecCopy(redsys->G, redsys->y));

    PetscCall(KSPSolve(kspA, redsys->F, tmp1));
    PetscCall(MatMult(BT, tmp1, tmp2));
    PetscCall(VecAXPY(redsys->y, -1, tmp2));

    // Solver for y = S^(-1)
    KSP ksp2;
	 PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp2));
    PetscCall(KSPSetOperators(ksp2, S, S));
    PetscCall(KSPSetType(ksp2, KSPMINRES));
    PetscCall(KSPSetInitialGuessNonzero(ksp2, PETSC_FALSE));
    PetscCall(KSPSetTolerances(ksp2, 1e-25, 1e-20, 10, 2000));
      
    PetscCall(KSPSolve(ksp2, redsys->y, redsys->y));

    // Compute x = A^{-1}(f-By)
    Vec tmp3;
    PetscCall(VecDuplicate(redsys->F, &redsys->x));
    PetscCall(VecCopy(redsys->F, redsys->x));

    PetscCall(VecDuplicate(redsys->F, &tmp3));

    PetscCall(MatMult(redsys->B, redsys->y, tmp3)); 
    PetscCall(VecAXPY(redsys->x, -1, tmp3));
    PetscCall(KSPSolve(kspA, redsys->x, redsys->x));

    return 0;
}
