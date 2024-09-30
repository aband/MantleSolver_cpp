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

    refArrayStokes_ = new int[br_->getDOF()];
    refArrayDarcy_  = new int[hdiv_->getDOF()];

    CreateRefMap(*br_, refArrayStokes_, mi, &bndryDOFStokes_);
    CreateRefMap(*hdiv_, refArrayDarcy_, mi, &bndryDOFDarcy_);

    Result_ = (ReducedSys *)malloc(sizeof(ReducedSys));

    return 1;
}

int Driver::SolveFlow(int maxIter, double tolUzawa){

    ParallelMatrixAssemble(mi, *basis_, myPhase, bndryStokesEssen_, reducedStokes_, 
                                                 bndryDarcyEssen_,  reducedDarcy_, 
                           &K_, *br_, *hdiv_ , mluseAdv_,

                           refArrayStokes_, refArrayDarcy_, bndryDOFStokes_, bndryDOFDarcy_);

    CreateLinearSys(reducedStokes_, M_*N_);
    CreateLinearSys(reducedDarcy_, M_*N_);

    CreateCoupledSystem(reducedStokes_, reducedDarcy_, Result_, &K_);

    CoupledUzawa(Result_, tolUzawa, maxIter);

    return 1;
}

int Driver::PrintFlowParallel(){

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

    CGNSPrepareParallel(&destStokes_sol, &destStokes_g, refArrayStokes_, mi,
                        ux, uy, *br_, *basis_);

    CGNSPrepareParallel(&destDarcy_sol, &destDarcy_g, refArrayDarcy_, mi,
                        vx, vy, *hdiv_, *basis_);

/*
    char stokesfile[] = "stokes.cgns";   
    CgnsArrayOutput(dmMesh,&globalmesh,ux,uy,mi.MPIlocalCellStart[0],
                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
                    mi.MPIlocalCellSize[1],stokesfile);

    char darcyfile[] = "darcy.cgns";    	
    CgnsArrayOutput(dmMesh,&globalmesh,vx,vy,mi.MPIlocalCellStart[0],
                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
                    mi.MPIlocalCellSize[1],darcyfile);
*/

    return 1;
}

int Driver::PrintCellCenterGrids(){

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

//        cout << HD << "  " << CD << "  " << global[1]  
//				 << global[1]*myPhase->pp->l0 << "   " 
//				 << lithoP << "  " << myPhase->pPtr->phi.mlt << endl;

        fprintf(fp, "%f ", myPhase->pPtr->phi.mlt);
    }
    fprintf(gridPorox, "\n");
    fprintf(gridPoroy, "\n");
	 fprintf(fp, "\n");}
    return 1;
}

inline bool exists_file (const std::string& name){
    struct stat buffer;

    return (stat (name.c_str(), &buffer) == 0);
}

int Driver::PrintPressure(){

    if (exists_file("gridCellX.dat") == 0){
        cout << "Yes " << endl;


    } 

    return 1;
}
