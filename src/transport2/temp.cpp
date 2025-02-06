#include "temp.h"

int printSol(int mark, Vec * global, const MeshInfo& mi){

    Vec temp = *global;

    const char * fieldname = "sol";

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    FILE * sol = fopen(filename,"w");

    for(int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for(int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double val;
        indice index {i,j};

        int nelem = FlatIndic(mi, index);

        PetscCall(VecGetValues(temp,1, &nelem, &val));

        fprintf(sol, "%e ", val);

    }fprintf(sol, "\n");}

    fclose(sol);

    return 1;
}

int printGrid(const MeshInfo& mi){

    FILE *gridPorox = fopen("gridCellX.dat", "w");
    FILE *gridPoroy = fopen("gridCellY.dat", "w");

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        vertex local {0.0,0.0};

        vector<vertex> corners = extractCorners(mi, {i,j});

        vertex global = GaussMapPointsFace(local, corners);

        fprintf(gridPorox,"%f ",global[0]);
        fprintf(gridPoroy,"%f ",global[1]);

    }
    fprintf(gridPorox, "\n");
    fprintf(gridPoroy, "\n");}

    fclose(gridPorox);
    fclose(gridPoroy);

    return 1;
}

std::string position(const indice& gcell){

    return "all";
}

int RK(double dt, int Nt, Vec * insol, const MeshInfo& mi, multilevel& ml,  mluse& use, DM dmu, DM dmmesh){

    Vec sol  = *insol;
    int event = 1;
    for (int t=0 ; t<Nt; t++){

        Vec flux;
        VecDuplicate(sol, &flux);

        getflux(mi, ml, use, &sol, &flux, dmu, dmmesh);

        VecAXPY(sol, -1*dt, flux);

        if (t%5 == 0){
        printSol(event,&sol,mi);
        event ++;
        }
    }
    printSol(event,&sol,mi);

    return 1;
}

int iRK(double dt, int Nt, Vec * insol, const MeshInfo& mi, multilevel& ml, mluse& use, DM dmu, DM dmmesh){

    int nelem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    Mat J;
    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           nelem, nelem, 
                           nelem, NULL, nelem, NULL, &J));
    PetscCall(MatSetUp((J)));

    KSP ksp;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetType(ksp, KSPGMRES));
    PetscCall(KSPSetTolerances(ksp, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE));

    Vec sol;
    VecDuplicate(*insol, &sol);
    VecCopy(*insol, sol);
    int event = 1;

    Vec flux;
    VecDuplicate(*insol, &flux);

    Vec previous;
    VecDuplicate(*insol, &previous);
    VecCopy(*insol, previous);


    Vec expsol;
    VecDuplicate(*insol, &expsol);
    VecCopy(*insol, expsol);

    Vec expflux;
    VecDuplicate(*insol, &expflux);

    for (int t=0; t<Nt; t++){

        double tol = 1.0;
        // Enter Newton's iteration
        int it = 0;
       
        //getflux(mi, ml, use, &sol, &flux, dmu, dmmesh);
        //VecAXPY(sol,-1.0*dt, flux);

        while (tol > 1e-7){
            Vec tmp1, tmp2;
            PetscCall(VecDuplicate(sol, &tmp1));
            PetscCall(VecDuplicate(sol, &tmp2));
            VecCopy(sol, tmp1);

            // 2. Compute function F(x) = x-previous + dt*f(x)
            // At the same time jacobian J(x) is computed
            // x stored in sol

            getall(mi, ml, use, &sol, &flux, &J, dmu, dmmesh, dt);

            VecAXPY(tmp1, -1.0, previous);
            VecAXPY(tmp1, 1.0, flux);

            // 3. Solve for J^-1(x)F(x)
            KSPSetOperators(ksp, J, J);

            KSPSolve(ksp, tmp1, tmp2);

            // 4. Update sol
            VecAXPY(sol, -1.0, tmp2);

            // 5. Compute 2nd norm of tmp2 and serve as tolerance indicator
            VecNorm(tmp2, NORM_2, &tol);

            cout << tol << endl;

            it ++;
        }

        cout << "Newton iteration count : " << it << endl;

        VecCopy(sol,previous);

        if (t%5 == 0){
        printSol(event,&sol,mi);
        event ++;
        }
    }

    printSol(event, &sol, mi);

    return 1;
}

int iRK2(double dt, int Nt, Vec * insol, MeshInfo* mi, multilevel* ml, mluse* use, DM dmu, DM dmmesh ){

    SNES snes;
    KSP  ksp;
    PC   pc;
    Vec  x, r; 
    Mat  J;

    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetType(snes, SNESNEWTONLS);

    VecDuplicate(*insol, &x);
    VecDuplicate(*insol, &r);

    VecCopy(*insol, x);

    int nelem = mi->MPIglobalCellSize[0] * mi->MPIglobalCellSize[1];

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           nelem, nelem, 
                           nelem, NULL, nelem, NULL, &J));
    PetscCall(MatSetUp(J));

    int event = 0;

    param * myparam = new param();
    myparam->mi = mi;
    myparam->ml = ml;
    myparam->use = use;
    myparam->dt = dt;
    myparam->dmu = dmu;
    myparam->dmmesh = dmmesh;
    myparam->previous = insol;

    Vec flux;
    VecDuplicate(*insol, &flux);
    myparam->flux = &flux;

    SNESSetFunction(snes, r, FormFunction, myparam);
    SNESSetJacobian(snes, J, J, FormJacobian, myparam);

    for (int t=0; t<Nt; t++){

        SNESSolve(snes, NULL, x);

        if (t%5 == 0){
        printSol(event,&x,*myparam->mi);
        event ++;
        }

        VecCopy(*myparam->previous, x);
    }

    return 1;
}

int getflux(const MeshInfo& mi, multilevel& ml, mluse& use, Vec * innow, Vec * influx, DM dmu, DM dmmesh){

    Vec now  = *innow; 
    Vec flux = *influx;

    Vec localu;

    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, now, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, now, INSERT_VALUES, localu); 

    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    double ** f;
    DMDAVecGetArray(dmu, flux, &f);

    // Update non linear weights with current cell-averaged solution
    ml.updatesigma(lu);

    Tensor<weights> allwgts;
    double h0 = sqrt((mi.L*mi.H)/(double)(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]));
    use.computeWgts(ml, mi, h0, allwgts);

    // Update edgeflux
    Tensor<double> horiedgeflux = Tensor<double>(2);
    horiedgeflux.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<double> vertedgeflux = Tensor<double>(2);
    vertedgeflux.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    updateEdgeFlux(vertedgeflux, horiedgeflux, mi, lu, use, ml, allwgts);

    // Loop through physical domain
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        f[j][i] = getcellflux(mi, {i,j}, vertedgeflux, horiedgeflux);
    }}

    DMDAVecRestoreArray(dmu, flux, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    return 1;
}

// Get jacobian as well
int getall(const MeshInfo& mi, multilevel& ml, mluse& use, 
           Vec * innow, Vec * influx, Mat *Jacobian, DM dmu, DM dmmesh, 
           const double& dt){

    Vec now  = *innow; 
    Vec flux = *influx;

    int nelem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           nelem, nelem, 
                           nelem, NULL, nelem, NULL, &(*Jacobian)));
    PetscCall(MatSetUp((*Jacobian)));

    Vec localu;

    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, now, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, now, INSERT_VALUES, localu); 

    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    double ** f;
    DMDAVecGetArray(dmu, flux, &f);

    double h0 = sqrt((mi.L*mi.H)/(double)(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]));

    // Update non linear weights with current cell-averaged solution
    ml.updateall(lu, h0, 1, 1e-4, mi);

    Tensor<weights> allwgts;
    use.computeWgts(ml, mi, h0, allwgts);

    // Update edgeflux
    Tensor<double> horiedgeflux = Tensor<double>(2);
    horiedgeflux.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<double> vertedgeflux = Tensor<double>(2);
    vertedgeflux.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    Tensor<derivative> horiedgefluxder = Tensor<derivative>(2);
    horiedgefluxder.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<derivative> vertedgefluxder = Tensor<derivative>(2);
    vertedgefluxder.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    updateEdgeFlux(vertedgeflux,    horiedgeflux, 
                   vertedgefluxder, horiedgefluxder,
                   mi, lu, use, ml, allwgts);

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double flux = 0.0;
        derivative dflux;

        getcellflux(mi, {i,j}, vertedgeflux, horiedgeflux, 
                    vertedgefluxder, horiedgefluxder, flux, dflux);

        f[j][i] = flux*dt;
        const int indexm = FlatIndic(mi, {i,j});

        for (const auto& it: dflux){
            const int indexn = it.first;
            const double val = it.second*dt;
            PetscCall(MatSetValues((*Jacobian), 1, &indexm, 1, &indexn, &val, ADD_VALUES));
        }

        const double val = 1.0;
        PetscCall(MatSetValues((*Jacobian), 1, &indexm, 1, &indexm, &val, ADD_VALUES));
 
    }}

    PetscCall(MatAssemblyBegin((*Jacobian), MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*Jacobian), MAT_FINAL_ASSEMBLY));
    DMDAVecRestoreArray(dmu, flux, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    return 1;
}


PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx){

    PetscFunctionBeginUser;
 
    param * user = (param *) ctx;

    //! Get local vector
    getflux((*user->mi), (*user->ml), (*user->use), &x, user->flux, user->dmu, user->dmmesh); 

    VecCopy(f, x);

    VecAXPY(f, -1.0, (*user->previous));
    VecAXPY(f, user->dt, (*user->flux));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void *ctx){

    PetscFunctionBeginUser;

    param * user = (param *)ctx;

    Vec flux;
    VecDuplicate(x, &flux);
    getall((*user->mi), (*user->ml), (*user->use), &x, &flux, &jac, user->dmu, user->dmmesh, user->dt); 

    PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
 
    PetscFunctionReturn(PETSC_SUCCESS);
}
