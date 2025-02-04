#include "print.h"

int printCellCenterGrid(const MeshInfo& mi){

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

int printCellAve(int mark, Vec * global, const MeshInfo& mi, const char * fieldname){

    Vec temp = *global;

    //const char * fieldname = "sol";

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

char * GetFilename(const char * fieldname, int mark){

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    return filename;
}

int Driver::PrintFlowEvent(int mark){

    // This function plots when distributed vectors have been scattered
    // This function should be called during the time stepping process
    // This function plot unscaled velocity of darcy velocity
    // This function also plots two effective velocity

    int nelemloc = mi.MPIlocalCellSize[0]*mi.MPIlocalCellSize[1];

    double * ux = (double *)malloc(sizeof(double)*nelemloc);
    double * uy = (double *)malloc(sizeof(double)*nelemloc);

    double * vx = (double *)malloc(sizeof(double)*nelemloc);
    double * vy = (double *)malloc(sizeof(double)*nelemloc);

    Vec stokesv;
    Vec darcyv;

    CreateScatterVec();

    CGNSPrepareParallel(&sresult_->vel_stokes, &sresult_->g_stokes, 
                        refArrayStokesEssen_, mi, ux, uy, *br_, *basis_, parameter); 

    CGNSPrepareParallel(&sresult_->vel_darcy, &sresult_->g_darcy,
                        refArrayDarcyEssen_, mi, vx, vy, *hdiv_, *basis_, parameter);

    quiverOutputEvent(ux,uy,vx,vy, mark, mi.MPIglobalCellSize[0], mi.MPIglobalCellSize[1]);

    return 1;
}

int quiverOutputEvent(double * ux, double * uy, double * vx, double * vy, int mark, int M, int N){

    // Separate velocity in x or y direction
    FILE * stokesVx = fopen(GetFilename("stokesVx",mark),"w");
    FILE * stokesVy = fopen(GetFilename("stokesVy",mark),"w");
    FILE * darcyVx  = fopen(GetFilename("darcyVx",mark),"w");
    FILE * darcyVy  = fopen(GetFilename("darcyVy",mark),"w");

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        fprintf(stokesVx, "%.21f ", ux[j*M+i]);
        fprintf(stokesVy, "%.21f ", uy[j*M+i]);
        fprintf(darcyVx, "%.21f ", vx[j*M+i]);
        fprintf(darcyVy, "%.21f ", vy[j*M+i]);
    }
    fprintf(stokesVx,"\n");
    fprintf(stokesVy,"\n");
    fprintf(darcyVx,"\n");
    fprintf(darcyVy,"\n");}

    fclose(stokesVx);
    fclose(stokesVy);
    fclose(darcyVx);
    fclose(darcyVy);

    return 1;
}

int Driver::PrintPhaseEvent(int mark){

    FILE *fp = fopen(GetFilename("porosity", mark), "w");
    FILE *ft = fopen(GetFilename("temperature", mark), "w");

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        vertex local {0.0,0.0};

        basis_->GetCorners(mi,{i,j});

        vertex global = GaussMapPointsFace(local, basis_->corners());
 
        double depth  = myPhase->pPtr->GetDepth(global[1], myPhase->pp->l0);  
        double lithoP = myPhase->pPtr->GetScaledLithoP(global[1]*(-1)*myPhase->pp->l0*0.6); 

        // Extract HD and CD from global solution vectors
        double HD, CD; // Get cell averaged values for approximation
        const int idx = FlatIndic(M, {i,j});
        PetscCall(VecGetValues(globalHD, 1, &idx, &HD));
        PetscCall(VecGetValues(globalCD, 1, &idx, &CD));

        myPhase->pPtr->evalPhase(HD, CD, lithoP);

        fprintf(fp, "%.21f ", myPhase->pPtr->phi.mlt);

        // Test ========================================================

        //double temp = AssignPorosity(global, myPhase->pp);
        //fprintf(fp, "%f ", temp);

        // =============================================================

        fprintf(ft, "%f ", myPhase->pPtr->TD);
    }
    fprintf(fp, "\n"); 
    fprintf(ft, "\n");}

    fclose(fp);
    fclose(ft);



    return 1;
}
