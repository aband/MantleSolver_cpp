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

int Driver::PrintFlowEventTransform(int mark){

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

    quiverOutputEventTransform(ux,uy,vx,vy, mark, mi.MPIglobalCellSize[0], mi.MPIglobalCellSize[1], V0);

    return 1;
}



int Driver::PrintEffVel(int mark, int side,
                        const Tensor<weights>& allwgtsHD, double ** lHD,
                        const Tensor<weights>& allwgtsCD, double ** lCD){

    // Print Effective computed at the gauss points of each edge
    // Only plot in selected direction
    FILE * effvx = fopen(GetFilename("effvelx", mark),"w");
    FILE * effvy = fopen(GetFilename("effvely", mark),"w");

    FILE * phasevx = fopen(GetFilename("phasevelx", mark),"w");
    FILE * phasevy = fopen(GetFilename("phasevely", mark),"w");

    FILE * solidvx = fopen(GetFilename("solidvelx", mark),"w");
    FILE * solidvy = fopen(GetFilename("solidvely", mark),"w");

    FILE * gaussgridx = fopen("gaussgridx.dat", "w");
    FILE * gaussgridy = fopen("gaussgridy.dat", "w");

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> gaussp;
    gaussp.resize(gpe.size());

    vertexSet edge;

    vector<vertex> effvel; effvel.resize(gaussp.size());
    vector<vertex> phasevel; phasevel.resize(gaussp.size());
    vector<vertex> solidvel; solidvel.resize(gaussp.size());
    vector<double> TDin; TDin.resize(gaussp.size());
    vector<double> TDout; TDout.resize(gaussp.size());
    vector<double> dTdHin; dTdHin.resize(gaussp.size());
    vector<double> dTdHout; dTdHout.resize(gaussp.size());
    vector<double> CDin; CDin.resize(gaussp.size());
    vector<double> CDout; CDout.resize(gaussp.size());
    vector<double> HDin; HDin.resize(gaussp.size());
    vector<double> HDout; HDout.resize(gaussp.size());
 
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice gcell {i,j};
        indice gcellout;
        vertexSet corners = extractCorners(mi, gcell);

        effvel.clear(); effvel.resize(gaussp.size());
        phasevel.clear(); phasevel.resize(gaussp.size());
        solidvel.clear(); solidvel.resize(gaussp.size());

        if (side==1){
        // vertical

            edge = {corners.at(3), corners.at(0)};

            for (int g=0; g<gpe.size(); g++){
                gaussp.at(g) = GaussMapPointsEdge({gpe[g]},edge);
            }   

            if (i==0){ // left side
                computeEffVel(gaussp, edge, gcell, allwgtsHD, lHD, allwgtsCD, lCD, effvel,phasevel,solidvel, TDin, dTdHin, CDin, HDin);
            } else {
                gcellout = gcell + mi.faceNormal[3];
                computeEffVel(gaussp, edge, gcell, gcellout, allwgtsHD, lHD, allwgtsCD, lCD, effvel,phasevel,solidvel, TDin, TDout, dTdHin, dTdHout, CDin, CDout, HDin, HDout);
            }

        } else if (side==2){
        // horizontal

            edge = {corners.at(0), corners.at(1)};

            for (int g=0; g<gpe.size(); g++){
                gaussp.at(g) = GaussMapPointsEdge({gpe[g]},edge);
            }   
            if (j==0){ // bottom side
                computeEffVel(gaussp, edge, gcell, allwgtsHD, lHD, allwgtsCD, lCD, effvel,phasevel,solidvel, TDin, dTdHin, CDin, HDin);
            } else {
                gcellout = gcell + mi.faceNormal[0];

                computeEffVel(gaussp, edge, gcell, gcellout, allwgtsHD, lHD, allwgtsCD, lCD, effvel,phasevel,solidvel, TDin, TDout, dTdHin, dTdHout, CDin, CDout, HDin, HDout);
            }

        } else {
            cout << "Pick a side. " << endl;
        }

        for (int g=0; g<gpe.size(); g++){

            fprintf(effvx, "%e ", effvel.at(g)[0]);
            fprintf(effvy, "%e ", effvel.at(g)[1]);

            fprintf(phasevx, "%e ", phasevel.at(g)[0]);
            fprintf(phasevy, "%e ", phasevel.at(g)[1]);

            fprintf(solidvx, "%e ", solidvel.at(g)[0]);
            fprintf(solidvy, "%e ", solidvel.at(g)[1]);

            fprintf(gaussgridx, "%e ", gaussp.at(g)[0]);
            fprintf(gaussgridy, "%e ", gaussp.at(g)[1]);
        }

    } fprintf(effvx, "\n ");
      fprintf(effvy, "\n ");
      fprintf(phasevx, "\n ");
      fprintf(phasevy, "\n ");
      fprintf(solidvx, "\n ");
      fprintf(solidvy, "\n ");

      fprintf(gaussgridx, "\n ");
      fprintf(gaussgridy, "\n ");}

    fclose(effvx);
    fclose(effvy);
    fclose(phasevx);
    fclose(phasevy);
    fclose(solidvx);
    fclose(solidvy);

    fclose(gaussgridx);
    fclose(gaussgridy);

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

int quiverOutputEventTransform(double * ux, double * uy, double * vx, double * vy, int mark, int M, int N, double V0){

    // Separate velocity in x or y direction
    FILE * stokesVx = fopen(GetFilename("stokesVx_transform",mark),"w");
    FILE * stokesVy = fopen(GetFilename("stokesVy_transform",mark),"w");
    FILE * darcyVx  = fopen(GetFilename("darcyVx_transform",mark),"w");
    FILE * darcyVy  = fopen(GetFilename("darcyVy_transform",mark),"w");

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        fprintf(stokesVx, "%.21f ", ux[j*M+i]*V0);
        fprintf(stokesVy, "%.21f ", uy[j*M+i]*V0);
        fprintf(darcyVx, "%.21f ", vx[j*M+i]*V0);
        fprintf(darcyVy, "%.21f ", vy[j*M+i]*V0);
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
    FILE *fphase = fopen(GetFilename("phase", mark), "w");
    FILE *ftm= fopen(GetFilename("meltT", mark), "w");
    FILE *fopx = fopen(GetFilename("opx", mark), "w");

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        vertex local {0.0,0.0};

        basis_->GetCorners(mi,{i,j});

        vertex global = GaussMapPointsFace(local, basis_->corners());
 
        double lithoP = myPhase->pPtr->GetStaticP(global[1]*(-1), myPhase->pp->l0); 

        // Extract HD and CD from global solution vectors
        double HD, CD; // Get cell averaged values for approximation
        const int idx = FlatIndic(M, {i,j});
        PetscCall(VecGetValues(globalHD, 1, &idx, &HD));
        PetscCall(VecGetValues(globalCD, 1, &idx, &CD));

        myPhase->pPtr->evalPhase(HD, CD, lithoP);

        fprintf(fp, "%e ", myPhase->pPtr->pc.phil);

        // Test ========================================================

        //double temp = AssignPorosity(global, myPhase->pp);
        //fprintf(fp, "%f ", temp);

        // =============================================================

        fprintf(ft, "%e ", myPhase->pPtr->pc.TDp);

        fprintf(fphase, "%d ", myPhase->pPtr->pc.region);

        fprintf(ftm, "%e ", myPhase->pPtr->GetTDp(myPhase->pPtr->TDe0, lithoP));

        fprintf(fopx, "%e ", myPhase->pPtr->pc.phi2);
    }
    fprintf(fp, "\n"); 
    fprintf(ft, "\n");
    fprintf(fphase, "\n");
	 fprintf(fopx, "\n");
    fprintf(ftm, "\n");}

    fclose(fp);
    fclose(ft);
    fclose(fphase);
    fclose(ftm);
    fclose(fopx);

    return 1;
}

int Driver::PrintPressureSerialApprox(int mark){
    // Print approximated pressure value in serial index system
    // This output is approximation of the pressure doing the following approximation:

    // 1. Pressure is represented by cell-averaged value
    // 2. Porosity computed using lithostatic pressure not actual pressure

    // Attention!! This approximated visual is for testing purpose only.

    // These two variables are solved in the reduced linear system
    Vec vectildeqf;
    Vec vecq;

    PetscCall(VecNestGetSubVec(Result_->y, 0, &vecq));   
    PetscCall(VecNestGetSubVec(Result_->y, 1, &vectildeqf));

    FILE * fqs = fopen(GetFilename("qs", mark), "w");
    FILE * fqf = fopen(GetFilename("qf", mark), "w");
    FILE * fps = fopen(GetFilename("ps", mark), "w");
    FILE * fpf = fopen(GetFilename("pf", mark), "w"); 

    FILE * fstokes = fopen(GetFilename("rawstokesq", mark), "w"); 
    FILE * fdarcy  = fopen(GetFilename("rawdarcyq", mark), "w");

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice globalcell {i,j};
        int nelem = FlatIndic(mi, globalcell);
        double q, tildeqf;

        PetscCall(VecGetValues(vectildeqf, 1, &nelem, &tildeqf));
        PetscCall(VecGetValues(vecq, 1, &nelem, &q));

        // Get Porosity
        vertex local {0.0,0.0};

        basis_->GetCorners(mi, globalcell);

        vertex global = GaussMapPointsFace(local, basis_->corners());
 
        double lithoP = myPhase->pPtr->GetStaticP(global[1]*(-1), myPhase->pp->l0); 

        // Extract HD and CD from global solution vectors
        double HD, CD; // Get cell averaged values for approximation
        const int idx = FlatIndic(mi, globalcell);
        PetscCall(VecGetValues(globalHD, 1, &idx, &HD));
        PetscCall(VecGetValues(globalCD, 1, &idx, &CD));

        myPhase->pPtr->evalPhase(HD, CD, lithoP);
        double phif = myPhase->pPtr->pc.phil;
        double coef = 0.0; 
        // Adjust phif
        if (phif > 2e-16) {
            coef = 1.0/sqrt(phif);
        }

        // Reterive original physical variables with physical units
        double qf = tildeqf * coef; 
        double qs = qf - 1.0/(1-phif)*(qf-q);
        double scale = myPhase->pp->rho_r * 10 * myPhase->pp->l0;

        qf *= scale;
        qs *= scale;

        double add = myPhase->pp->rho_f * 10 * global[1]*myPhase->pp->l0;

        double pf = qf + add;
        double ps = qs + add; 

        fprintf(fqs, "%e ", qs);
        fprintf(fqf, "%e ", qf);
        fprintf(fps, "%e ", ps);
        fprintf(fpf, "%e ", pf);
        fprintf(fstokes, "%e ", q);
        fprintf(fdarcy, "%e ", tildeqf);
    }fprintf(fqs, "\n");
     fprintf(fps, "\n");
     fprintf(fqf, "\n");
     fprintf(fps, "\n");
     fprintf(fstokes, "\n");
     fprintf(fdarcy, "\n");}

    return 1;
}
