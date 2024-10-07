#include "driver.h"

int Driver::PrepareFlow(){

    basis_ = new basis();
    hdiv_  = new Hdivmixed();
    br_    = new BRMixed();

    br_->ComputeTotalDOF(mi);
    hdiv_->ComputeTotalDOF(mi);

    MarkBndryDOFStokes(bndryStokesEssen_, bndryStokesNatur_, mi, *basis_, *br_, myPhase->pp);
    MarkBndryDOFDarcy(bndryDarcyEssen_, bndryDarcyNatur_, mi, *basis_, *hdiv_, myPhase->pp);

    reducedDarcy_ = (ReducedSys *)malloc(sizeof(ReducedSys));
    reducedStokes_ = (ReducedSys *)malloc(sizeof(ReducedSys));

    refArrayStokesEssen_ = new int[br_->getDOF()];
    refArrayDarcyEssen_  = new int[hdiv_->getDOF()];

//    CreateRefMap(*br_, refArrayStokesEssen_, mi, &bndryDOFStokes_);
//    CreateRefMap(*hdiv_, refArrayDarcyEssen_, mi, &bndryDOFDarcy_);

    CreateRefMap(*br_  , mi, refArrayStokesEssen_, refArrayStokesNatur_, &bndryDOFStokes_, &bndryDOFStokesNatur_);
    CreateRefMap(*hdiv_, mi, refArrayDarcyEssen_ , refArrayDarcyNatur_ , &bndryDOFDarcy_ , &bndryDOFDarcyNatur_ );

    Result_ = (ReducedSys *)malloc(sizeof(ReducedSys));

    return 1;
}

int Driver::SolveFlow(int maxIter, double tolUzawa){

    ParallelMatrixAssemble(mi, *basis_, myPhase, bndryStokesEssen_, reducedStokes_, 
                                                 bndryDarcyEssen_,  reducedDarcy_, 
                           &K_, *br_, *hdiv_ , mluseAdv_,

                           refArrayStokesEssen_, refArrayDarcyEssen_, bndryDOFStokes_, bndryDOFDarcy_);

    CreateLinearSys(reducedStokes_, M_*N_);
    CreateLinearSys(reducedDarcy_, M_*N_);

    CreateCoupledSystem(reducedStokes_, reducedDarcy_, Result_, &K_);

    CoupledUzawa(Result_, tolUzawa, maxIter);

    return 1;
}

inline int quiverOutputSerial(double * ux, double * uy, double *vx, double *vy, int M, int N){

    // Output of velocity on vertex without using cgns format
    // Working in serial case only

    FILE * stokesVx = fopen("stokesVx.dat", "w");
    FILE * stokesVy = fopen("stokesVy.dat", "w");
    FILE * darcyVx  = fopen("darcyVx.dat", "w");
    FILE * darcyVy  = fopen("darcyVy.dat", "w"); 

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        fprintf(stokesVx, "%f ", ux[j*M+i]);
        fprintf(stokesVy, "%f ", uy[j*M+i]);
        fprintf(darcyVx, "%f ", vx[j*M+i]);
        fprintf(darcyVy, "%f ", vy[j*M+i]);
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

int Driver::PrintFlow(){

    // Reconstruct values at the centroid of cells
    int nelemloc = mi.MPIlocalCellSize[0]*mi.MPIlocalCellSize[1];
    double * ux = (double *)malloc(sizeof(double)*nelemloc);
    double * uy = (double *)malloc(sizeof(double)*nelemloc);

    double * vx = (double *)malloc(sizeof(double)*nelemloc);
    double * vy = (double *)malloc(sizeof(double)*nelemloc);
 
    Vec stokesv;
    Vec darcyv;

    PetscCall(VecNestGetSubVec(Result_->x, 0, &stokesv));
    PetscCall(VecNestGetSubVec(Result_->x, 1, &darcyv));

    Vec destStokes_sol, destStokes_g;
    Vec destDarcy_sol, destDarcy_g;

    SolScatAll(&stokesv, &reducedStokes_->g, 
               &destStokes_sol, &destStokes_g);  

    SolScatAll(&darcyv, &reducedDarcy_->g, 
               &destDarcy_sol, &destDarcy_g);  

    CGNSPrepareParallel(&destStokes_sol, &destStokes_g, refArrayStokesEssen_, mi,
                        ux, uy, *br_, *basis_);

    CGNSPrepareParallel(&destDarcy_sol, &destDarcy_g, refArrayDarcyEssen_, mi,
                        vx, vy, *hdiv_, *basis_);

#ifdef CGNS_OUT
    char stokesfile[] = "stokes.cgns";   
    CgnsArrayOutput(dmMesh,&globalmesh,ux,uy,mi.MPIlocalCellStart[0],
                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
                    mi.MPIlocalCellSize[1],stokesfile);

    char darcyfile[] = "darcy.cgns";    	
    CgnsArrayOutput(dmMesh,&globalmesh,vx,vy,mi.MPIlocalCellStart[0],
                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
                    mi.MPIlocalCellSize[1],darcyfile);
#endif

#ifndef CGNS_OUT
    quiverOutputSerial(ux,uy,vx,vy,M_,N_);
#endif

    return 1;
}

int Driver::PrintPorosity(){

    // Cell centered grid
    FILE *gridPorox = fopen("gridCellX.dat", "w");
    FILE *gridPoroy = fopen("gridCellY.dat", "w");

    FILE *fp = fopen("porosity.dat","w");

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        vertex local {0.0,0.0};

        basis_->GetCorners(mi,{i,j});

        vertex global = GaussMapPointsFace(local, basis_->corners());

        fprintf(gridPorox,"%f ",global[0]);
        fprintf(gridPoroy,"%f ",global[1]);
 
        double HD = InitHD(global,{myPhase->pp->l0,0.0});
        double CD = InitCD(global,{0.0,0.0});

        double lithoP = myPhase->pPtr->GetScaledLithoP(global[1]*(-1)*myPhase->pp->l0*0.6); 

        myPhase->pPtr->evalPhase(HD, CD, lithoP);

        fprintf(fp, "%f ", myPhase->pPtr->phi.mlt);
    }
    fprintf(gridPorox, "\n");
    fprintf(gridPoroy, "\n");
    fprintf(fp, "\n");}

    fclose(gridPorox);
    fclose(gridPoroy);
    fclose(fp);

    return 1;
}

inline bool exists_file (const std::string& name){
    struct stat buffer;

    return (stat (name.c_str(), &buffer) == 0);
}

int Driver::PrintPressureConstant(){

    // Print piece wise constant pressure
    if (exists_file("gridCellX.dat") == 0 || exists_file("gridCellY.dat") == 0){
        // Checking for cell centered grid file
        // These two grid files are output along with printflow 
        cout << "No Grid File " << endl;
    } 
	 
    FILE * darcyp  = fopen("darcyp.dat","w");
    FILE * stokesp = fopen("stokesp.dat","w");

    Vec stokesP;
    Vec darcyP;

    PetscCall(VecNestGetSubVec(Result_->y, 0, &stokesP));
    PetscCall(VecNestGetSubVec(Result_->y, 1, &darcyP));

    int istart = mi.MPIlocalCellStart[0];
    int jstart = mi.MPIlocalCellStart[1];

    for (int j=jstart; j<jstart + mi.MPIlocalCellSize[1]; j++){
    for (int i=istart; i<istart + mi.MPIlocalCellSize[0]; i++){
        indice global {i,j};
        int nelem = FlatIndic(mi, global);
        // Extract pressure from Vec
        double dp, sp;

        PetscCall(VecGetValues(stokesP,1, &nelem, &sp));
        PetscCall(VecGetValues(darcyP,1, &nelem, &dp));

        fprintf(darcyp, "%f ", dp);
        fprintf(stokesp, "%f ", sp);

    }fprintf(darcyp, "\n");
     fprintf(stokesp, "\n");}

    fclose(stokesp);
    fclose(darcyp);

    return 1;
}
