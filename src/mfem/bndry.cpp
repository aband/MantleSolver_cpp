PetscErrorCode AssignValuesRHS(int NS, int ND, int Nelem,
                               Vec * A, Vec * B,
                               RHSVector * rhsv,
                               const unordered_map<int,int>& bndryStokes,
                               const unordered_map<int,int>& bndryDarcy){

    PetscFunctionBeginUser;

    Vec a, b;

    a = *A;
    b = *B;

    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, ND, &rhsv->ad));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, NS, &rhsv->bs));

    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, Nelem, &rhsv->qd));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, Nelem, &rhsv->qs));

    // Assign local arrays
    double *localad;
    double *localbs;
    double *localqd;
    double *localqs;

    double *localsource;

    PetscCall(VecGetArray(rhsv->ad, &localad));
    PetscCall(VecGetArray(rhsv->bs, &localbs));
    PetscCall(VecGetArray(rhsv->qd, &localqd));
    PetscCall(VecGetArray(rhsv->qs, &localqs));

    PetscCall(VecGetArray(rhsv->source, &localsource));

    // Initialize local arrays
    for (int k=0; k<ND; k++){localad[k] = 0.0;}
    for (int k=0; k<NS; k++){localbs[k] = 0.0;}
    for (int k=0; k<Nelem; k++) {localqd[k] = 0.0; localqs[k] = 0.0;}




    PetscFunctionReturn(0);
}

int ComputeBndryValsStokes(const unordered_map<int, int>& bndryStokes){




    return 0;
}

double * ComputeBndryValsDarcy(const unordered_map<int, int>& bndryDarcy){

}

PetscErrorCode CreateRHS(const MeshInfo& mi,
                         basis& basis_,
                         Hdivmixed& hdiv_,
                         BRMixed& br_,
                         PhysProperty * physproperty,
                         RHSVector * rhsv){

    // Create right hand side vector
    PetscFunctionBeginUser;

    // Create auxilliary vectors
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, ));


    PetscFunctionReturn(0);
}
