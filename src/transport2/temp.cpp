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

// Initial consition
double func(const vertex& point,
            const vector<double>& param){

    // Initial condition

    // Initialize with simple Reimann shock and rarefaction function
    // time inputed as param[0] 

    if (point[0] < 0.5 || point[0] >=(0.5*param[0]+1.5)){
        return 0;
    } else if (point[0]>=0.5 && point[0]<param[0]+0.5){
        return (point[0]-0.5)/param[0];
    } else if (point[0]>=point[0]+0.5 || point[0] <0.5*param[0]+1.5){
        return 1;
    } else {
        return 0;
    }

}

// Burgers for testing
double advfunc(const double& u, 
               const vertex& vel, const vertex& unitnormal){

    // A Burgers type flux

    return u*u/2.0 *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
}

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work){

    // compute df/du = df/dR * dR/du

    double direction = u*(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);

    unordered_map_arithmetic(work, du, std::plus<double>(), direction, std::multiplies<double>());

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

    KSP ksp;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetType(ksp, KSPGMRES));
    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE));

    Mat J;

    Vec sol = *insol;
    int event = 1;

    Vec flux;
    VecDuplicate(sol, &flux);

    Vec previous;
    VecDuplicate(sol, &previous);
    VecCopy(sol, previous);

    double tol = 0.0;

    for (int t=0; t<Nt; t++){

        // Newton's iteration 
        // 1. Get initial guess x0 = u-dtF
        // store x0 in sol
        getflux(mi, ml, use, &sol, &flux, dmu, dmmesh);
        VecAXPY(sol, -1*dt, flux);

        // Enter Newton's iteration
        while (tol > 1e-10){

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
            VecNorm(sol, NORM_2, &tol);
        }

        VecCopy(sol,previous);

        if (t%5 == 0){
        printSol(event,&sol,mi);
        event ++;
        }
    }

    printSol(event, &sol, mi);

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

    Mat J = *Jacobian;

    int nelem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           nelem, nelem, 
                           nelem, NULL, nelem, NULL, &J));
    PetscCall(MatSetUp(J));

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

        const int indexn = FlatIndic(mi, {i,j});

        for (const auto& it: dflux){
            const int indexm = it.first;
            const double val = it.second*dt;

            PetscCall(MatSetValues(J, 1, &indexm, 1, &indexn, &val, ADD_VALUES));
        }

        const double val = 1.0;
        PetscCall(MatSetValues(J, 1, &indexn, 1, &indexn, &val, ADD_VALUES));
 
    }}

    PetscCall(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));

    DMDAVecRestoreArray(dmu, flux, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    return 1;
}
