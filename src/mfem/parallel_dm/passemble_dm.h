#ifndef PASSEMBLE_DM_H_
#define PASSEMBLE_DM_H_

#include "passemble.h"

inline int PrepareFullSys(FullSys * sys, DM dmv, DM dmp, int Bdnx, int Bonz){

    // Two different kinds of data management
    // dmv used for velocity space, a DMStag object
    // dmp used for pressure space, a DMDA object

    PetscCall(DMCreateMatrix(dmv, &sys->A));

    int sizeA;
    PetscCall(MatGetSizes(sys->A, NULL, NULL, sizeA, NULL));

    int sizeC;
    PetscCall(MatGetSizes(sys->C, NULL, NULL, sizeC, NULL));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           sizeA, sizeC, Bdnz, NULL, Bonz, NULL, &sys->B));
    PetscCall(MatSetUp(sys->B));

    return 0;
}

#endif
