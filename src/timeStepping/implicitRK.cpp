// Parallel implicit Runge Kutta based upon Petsc

typedef struct{
    DM dm;
    vector<WenoReconstruction *> wr;
    MeshInfo mi;
    int stencil_count;
} Ctx;

PetscErrorCode FormJacobianSNES(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void * ctx){
    PetscErrorCode    ierr;
    PetscFunctionBeginUser;




    PetscFunctionReturn(0);
}

PetscErrorCode ImplicitRungeKutta(const vector<double>& c, 
                                  const vector<double>& bT, 
                                  const vector<double>& A,
                                  int totalStage;
                                  double h, double T, 
                                  const MeshInfo& mi, ... int matSize){

    PetscErrorCode    ierr;
    PetscFunctionBeginUser;            

    assert(std::accumulate(bT.begin(), bT.end(), decltype(bT)::value_type(0)) == 1.0);

    // Solving the first linear system z - kron(A,I) F(z) = 0 first.

    SNES snes;
    Ctx  ctx; 

    // Set up ctx data
    ctx.wr = wr;
    ctx.mi = mi;
    ctx.stencil_count = stencil_count;

    Mat iRK;
    ierr = MatCreateBAIJ(PETSC_COMM_WORLD, size, PETSC_DECIDE, PETSC_DECIDE, 
                         totalStage*matSize, totalStage*matSize, 0, NULL, 1, NULL, &iRK);CHKERRQ(ierr);

    ierr = SNESSetJacobian(snes, iRK, iRK, FormJacobianSNES, &ctx);CHKERRQ(ierr);


    PetscFunctionReturn(0);
}
