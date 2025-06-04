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

    // Schur Complement (C - BTA-1B)
    Mat S;
    MatCreateSchurComplement(redsys->M, redsys->M, redsys->B, BT, redsys->C, &S);
    MatSchurComplementGetKSP(S, &kspA);
    
    KSP kspSchur;
    KSPSetOperators(kspSchur, S, S);

    // Compute g-BTA-1f
    Vec tmp1, tmp2, tmp3;
    PetscCall(VecDuplicate(redsys->F, &tmp1));
    PetscCall(VecDuplicate(redsys->F, &tmp2));
    PetscCall(VecDuplicate(redsys->F, &tmp3));
 
    PetscCall(KSPSolve(kspA, redsys->F, tmp1));
    PetscCall(MatMult(BT, tmp1, tmp2));
    PetscCall(VecAXPY());

    return 0;
}
