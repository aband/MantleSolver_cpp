PetscErrorCode ParallelMatrixAssemble(DM dmstag){

    PetscInt startx, starty, nx, ny, M, N;

    PetscFunctionBeginUser;

    PetscCall(DMStagGetCorners(dmstag, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
    PetscCall(DMStagGetGlobalSizes(dmstag, &M, &N, NULL););

    for (int j=starty; j<starty+ny; j++){
    for (int i=startx; i<startx+nx; i++){


    }}

    return PETSC_SUCCESS;
}
