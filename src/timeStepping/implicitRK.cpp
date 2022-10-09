// Parallel implicit Runge Kutta based upon Petsc
#include "../../include/implicitRK.h"

PetscErrorCode DrawMat(Mat V, char * myfile)
{
    PetscFunctionBeginUser;

    FILE *f = fopen(myfile,"w");

    if (f == NULL){
        printf("Error opening file !\n");
        exit(1);
    }

    int mm,nn;
    MatGetSize(V,&nn,&mm);

    for (int j=nn-1; j>0; j--){
    for (int i=0; i<mm; i++){
        double a;
        MatGetValues(V,1,&j,1,&i,&a);
        fprintf(f,"%f ",a);
    }fprintf(f,"\n ");}

    fclose(f);

    PetscFunctionReturn(0);
}

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

PetscErrorCode FormJacobianIEULER(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void * ctx){

    PetscErrorCode    ierr;
    Ctx *user = (Ctx*)ctx;
    DM  dm = (DM)user->dm;
    PetscInt M,N,xs,ys,xm,ym,stencilwidth;
    PetscFunctionBeginUser;

/*
 * Set up necessary variables for computation of jacobian
 */
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

    user->mi.localval = lu;

    // Update corresponding non linear weights
    for (int s=0; s<user->stencil_count; s++){
        wr[s]->ComputeNonlinWeights(user->mi);
    }

    int rstart, rend;

    MatGetOwnershipRange(J, &rstart, &rend);

    // assume periodic boundary condition
    for (int row = rstart; row<rend; row++){

        int nWr = (row/M+1)*(M+2) + (row%M)+1;

        index_set gIndex = wr[nWr]->GetGlobalCellIndexStencil(user->mi);

        int vertxx = row%M;
        int vertxy = row/M;

        point_index target {vertxx+user->mi.ghost_vertx[0], 
                            vertxy+user->mi.ghost_vertx[1]};

        vector<double> deriv = DerivLaxFriedrichFlux(user->mi, time, target, wr, funcX, funcY, dfuncX, dfuncY); 

        for (int d=0; d<deriv.size(); d++){
            // Define placement of elements
            // Periodic boundary for now

            PetscInt col = (gIndex[d][1]*M + gIndex[d][0] + M*N)%(M*N);

            PetscScalar val = -1.0*deriv[d];

            //if (col  == row){
            //    val = 1.0 - deriv[d];
            //} else {
            //    val = -1.0 * deriv[d];
            //}

            ierr = MatSetValue(J,row,col,val,ADD_VALUES) ; CHKERRQ(ierr);

        }

    }

    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (J != Jp){
        ierr = MatAssemblyBegin(Jp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Jp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

    ierr = DMDAVecRestoreArray(dm, localu, &lu);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &localu);CHKERRQ(ierr);

    char file[] = "matrix.data";
    DrawMat(J, file);

    PetscFunctionReturn(0);
}

PetscErrorCode FormJacobianIRK(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void * ctx){

    PetscErrorCode    ierr;
    PetscFunctionBeginUser;




    PetscFunctionReturn(0);
}

PetscErrorCode MPIImplicitRungeKutta(const vector<double>& c, 
                                  const vector<double>& bT, 
                                  const vector<double>& A,
                                  int totalStage,
                                  double h, double T, 
                                  void * ctx){

    PetscErrorCode    ierr;
    PetscMPIInt       size, rank;
    PetscFunctionBeginUser;            

    //assert(std::accumulate(bT.begin(), bT.end(), decltype(bT)::value_type(0)) == 1.0);
    assert(A.size() == totalStage*totalStage);

    // Convert vector to array
    // Check if A is invertable or not
    int N = totalStage;
    double *arrayA = new double[N];
    int *IPIV = new int[N];
    int LWORK = N*N;
    double *inverseA = new double[LWORK];

    copy(A.begin(), A.end(), arrayA); 
    //dgetri(&N, arrayA, &N, IPIV, inverseA, &LWORK, &ierr);
    assert(ierr == 0);

    // Get MPI Info
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Solving the first linear system z - kron(A,I) F(z) = 0 first.
    SNES snes;

    //int totalSize = totalStage*matSize;

    Mat iRK;
    //ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
    //                    totalSize, totalSize, , NULL, 1, NULL, &iRK);CHKERRQ(ierr);


    //ierr = SNESSetJacobian(snes, iRK, iRK, FormJacobianSNES, &ctx);CHKERRQ(ierr);
    
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

    int localSize = PetscFloorReal(fullSize / size); 

    if (size == rank){
        localSize = fullSize - localSize*(size-1);
    }

    // MPI size should be larger than max polynomial order

    //ierr = VecCreateMPIWithArray(); CHKERRQ(ierr);

    //ierr = MatCreateAIJ(PETSC_COMM_WORLD, localsize, ); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}

PetscErrorCode MPIImplicitEuler(double h, double T, void * ctx){

    PetscFunctionBeginUser;

    PetscFunctionReturn(0);
}

/*
 * A sequential function used to calculate implicit euler
 * time stepping for reference. Implicit Euler can be done
 * with the help of time stepping onject TS.
 */

PetscErrorCode SeqImplicitEuler(int stencil_count, vector<double>& linWeights, vector<int *>& rangex, vector<int *>& rangey, MeshInfo& mi, DM dmu, double T, double dt, Vec globalu){

    PetscErrorCode    ierr;
    PetscMPIInt       size, rank;
    PetscFunctionBeginUser;
   
    // Get MPI Info
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    assert(size == 1); // This should be a sequential code

    PetscInt M,N,xs,ys,xm,ym,stencilwidth;
    ierr = DMDAGetInfo(dmu, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Create weno reconstruction class
    vector<WenoReconstruction *> wr;
    wr.resize(stencil_count); // Should be local stencil count instead of global count
    
    for (int s=0; s<stencil_count; s++){
        int shiftj = s/(M+2)-1;
        int shifti = s%(M+2)-1;
        valarray<int> target = {shifti, shiftj};
        wr[s] = new WenoReconstruction(mi,linWeights,rangex,rangey,target);
    }

    // Set up time stepping
    TS   ts;
    SNES snes;
    Ctx  ctx;

    // Time stepping with TS object
    ctx.wr = wr;
    ctx.dm = dmu;
    ctx.mi = mi;
    ctx.stencil_count = stencil_count;

    // Allocate space for Jacobian computation
    Mat A;
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, M*N, M*N);
    MatSetUp(A);

    // Set up snes
    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetType(snes, SNESNGMRES);

	 TSCreate(PETSC_COMM_WORLD, &ts);
	 TSSetProblemType(ts,TS_NONLINEAR);
	 TSSetType(ts, TSBEULER);

	 TSSetMaxTime(ts,T);
	 TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);
	 TSSetDM(ts,dmu);

	 // Customize nonlinear lu[j][i]
	 TSGetSNES(ts,&snes);
	 TSSetTimeStep(ts,dt);
	 TSSetSolution(ts,globalu);

	 TSSetRHSFunction(ts, globalu, FormFunction, &ctx);

    TSSetRHSJacobian(ts, A, A, FormJacobianIEULER, &ctx);

	 cout << "Time stepping started." << endl;
	 cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

	 TSSolve(ts,globalu);

	 cout << "Time stepping ended." << endl;

	 // ==========================================================================

    PetscFunctionReturn(0);
}
