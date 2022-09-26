// Parallel implicit Runge Kutta based upon Petsc

typedef struct{
    DM dm;
    vector<WenoReconstruction *> wr;
    MeshInfo mi;
    int stencil_count;
} Ctx;

// Form function using pure transport flux
// the pointer to totalflux function should be
// reconsidered
PetscErrorCode FormFunction(TS ts, PetscReal time, Vec U, Vec F, void * ctx){

    PetscErrorCode ierr;
    Ctx *user = (Ctx*)ctx;
    DM  dm = (DM)user->dm;
    PetscInt M,N,xs,ys,xm,ym,stencilwidth;
    PetscFunctionBeginUser;

    vector<WenoReconstruction*>& wr = user->wr;

    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                                CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Get local vector
    Vec localu;
    DMGetLocalVector(dm, &localu);

    DMGlobalToLocalBegin(dm, U, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dm, U, INSERT_VALUES, localu);

    // It can be changed later to not be double
    double  ** lu;
    DMDAVecGetArray(dm, localu, &lu);

    double ** f;
    DMDAVecGetArray(dm, F, &f);

    user->mi.localval = lu;

    // Calculate corresponding non linear weights
    for (int s=0; s<user->stencil_count; s++){
        wr[s]->ComputeNonlinWeights(user->mi);
    }

    int offset = 4;

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+ym; i++){
        if (j<offset || i<offset || j>N-offset || i>M-offset){
            f[j][i] = lu[j][i];
        } else {
            point_index target {i-xs+user->mi.ghost_vertx[0], j-ys+user->mi.ghost_vertx[1]};
            double temp = 0.0;
            for (int pos = 0; pos<4; pos++){
                temp -= 1.0/(wr[(j-ys+1)*(xm+2)+(i-xs+1)]->Geth()*wr[(j-ys+1)*(xm+2)+(i-xs+1)]->Geth())
                        * TotalFlux(user->mi, pos, time, target, wr, funcX, funcY, dfuncX, dfuncY);
            }
            f[j][i] = temp;
        }
    }}

    DMDAVecRestoreArray(dm, F, &f);
    DMDAVecRestoreArray(dm, localu, &lu);
    DMRestoreLocalVector(dm, &localu);

    PetscFunctionReturn(0);
}

PetscErrorCode FormJacobianIRK(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void * ctx){

    PetscErrorCode    ierr;
    PetscFunctionBeginUser;




    PetscFunctionReturn(0);
}

PetscErrorCode FormJacobianIEULER(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void * ctx){

    PetscErrorCode    ierr;
    PetscFunctionBeginUser;

/*
 *
 * Set up necessary variables for computation of jacobian
 *
 */
    int M = size of U;  

    int layer = length of ghost layer;

    // assume periodic boundary condition
    for (int i=0; i<M; i++){
        if (i<layer || i>M-layer){
            // The situation where periodic boundary condition kicks in

        } else {
            vector<double> 
        }
    }


    PetscFunctionReturn(0);
}

PetscErrorCode MPIImplicitRungeKutta(const vector<double>& c, 
                                  const vector<double>& bT, 
                                  const vector<double>& A,
                                  int totalStage;
                                  double h, double T, 
                                  void * ctx){

    PetscErrorCode    ierr;
    PetscMPIInt       size, rank;
    PetscFunctionBeginUser;            

    assert(std::accumulate(bT.begin(), bT.end(), decltype(bT)::value_type(0)) == 1.0);
    assert(A.size() == totalStage*totalStage);

    // Convert vector to array
    // Check if A is invertable or not
    int N = totalStage;
    double *arrayA = new double[N];
    int *IPIV = new int[N];
    int LWORK = N*N;
    double *inverseA = new double[LWORK];

    copy(A.begin(), A.end(), arrayA); 
    dgetri(&N, arrayA, &N, IPIV, inverseA, &LWORK, &ierr);
    assert(ierr == 0);

    // Get MPI Info
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Solving the first linear system z - kron(A,I) F(z) = 0 first.
    SNES snes;

    int totalSize = totalStage*matSize;

    Mat iRK;
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                        totalSize, totalSize, , NULL, 1, NULL, &iRK);CHKERRQ(ierr);


    ierr = SNESSetJacobian(snes, iRK, iRK, FormJacobianSNES, &ctx);CHKERRQ(ierr);
    
    // Delete declared dynamic memory array
    // Be careful with deleting array
    delete[] IPIV;
    delete[] inverseA;
    delete[] arrayA; 
    PetscFunctionReturn(0);
}

// Simplist implicit time stepping for checking the correctness of Jacobian
PetscErrorCode MPIPseudoImplicitEuler(int fullSize, double h, double T, void * ctx){

    // Using Pseudo Jacobian instead of full Jacobian
    PetscErrorCode    ierr;
    PetscMPIInt       size, rank;
    PetscFunctionBeginUser;
   
    // Get MPI Info
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


    SNES snes;

    // Create Jacobian Mat
    Mat iEuler;
    Vec U;

    int localSize = PetscFloorReal(FullSize / size); 

    if (size == rank){
        localSize = fullSize - localSize*(size-1);
    }

    // MPI size should be larger than max polynomial order

    ierr = VecCreateMPIWithArray(); CHKRRQ(ierr);

    ierr = MatCreateAIJ(PETSC_COMM_WORLD, localsize, ); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}

PetscErrorCode MPIImplicitEuler(double h, double T, void * ctx){

    PetscFunctionBeginUser;

    PetscFunctionReturn(0);
}

// A sequential function used to calculate implicit euler
// time stepping for reference
PetscErrorCode SeqImplicitEuler(){

    PetscErrorCode    ierr;
    PetscMPIInt       size, rank;
    PetscFunctionBeginUser;
   
    // Get MPI Info
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    assert(size == 1); // This should be a sequential code

    

    PetscFunctionReturn(0);
}

