#include "driver.h"

int Driver::PrepareFlow(){

    basis_ = new basis();
    hdiv_  = new Hdivmixed();
    br_    = new BRMixed();

    br_->ComputeTotalDOF(mi);
    hdiv_->ComputeTotalDOF(mi);
  
    // Define boundary parameter
    parameter.push_back(-0.2);

    MarkBndryDOFStokes(bndryStokesEssen_, bndryStokesNatur_, mi, *basis_, *br_, myPhase->pp, parameter);
    MarkBndryDOFDarcy(bndryDarcyEssen_, bndryDarcyNatur_, mi, *basis_, *hdiv_, myPhase->pp);

    reducedDarcy_ = (ReducedSys *)malloc(sizeof(ReducedSys));
    reducedStokes_ = (ReducedSys *)malloc(sizeof(ReducedSys));

    refArrayStokesEssen_ = new int[br_->getDOF()];
    refArrayDarcyEssen_  = new int[hdiv_->getDOF()];

//    CreateRefMap(*br_, refArrayStokesEssen_, mi, &bndryDOFStokes_);
//    CreateRefMap(*hdiv_, refArrayDarcyEssen_, mi, &bndryDOFDarcy_);

    CreateRefMap(*br_  , mi, refArrayStokesEssen_, refArrayStokesNatur_, &bndryDOFStokes_, &bndryDOFStokesNatur_, parameter);
    CreateRefMap(*hdiv_, mi, refArrayDarcyEssen_ , refArrayDarcyNatur_ , &bndryDOFDarcy_ , &bndryDOFDarcyNatur_ , parameter);

    Result_ = (ReducedSys *)malloc(sizeof(ReducedSys));

    sresult_ = (ScatterResult *)malloc(sizeof(ScatterResult));

    return 1;
}

int Driver::SolveFlow(int maxIter, double tolUzawa){

    ParallelMatrixAssemble(mi, *basis_, myPhase, bndryStokesEssen_, reducedStokes_, 
                                                 bndryDarcyEssen_,  reducedDarcy_, 
                           &K_, *br_, *hdiv_ , mluseAdv_,

                           refArrayStokesEssen_, refArrayDarcyEssen_, bndryDOFStokes_, bndryDOFDarcy_, parameter);

    CreateLinearSys(reducedStokes_, M_*N_);
    CreateLinearSys(reducedDarcy_, M_*N_);

    CreateCoupledSystem(reducedStokes_, reducedDarcy_, Result_, &K_);

    CoupledUzawa(Result_, tolUzawa, maxIter);

    return 1;
}

int Driver::CreateScatterVec(){

    // Scatter distributed vector to all processors
    Vec stokesv, darcyv;
    PetscCall(VecNestGetSubVec(Result_->x, 0, &stokesv));
    PetscCall(VecNestGetSubVec(Result_->x, 1, &darcyv));

    SolScatAll(&stokesv, &reducedStokes_->g, 
               &sresult_->vel_stokes, &sresult_->g_stokes);  

    SolScatAll(&darcyv, &reducedDarcy_->g, 
               &sresult_->vel_darcy, &sresult_->g_darcy);  

    return 1;
}

int Driver::quiverOutputEvent(double * ux, double * uy, double * vx, double * vy){

    // Separate velocity in x or y direction
    FILE * stokesVx = fopen(GetFilename("stokesVx"),"w");
    FILE * stokesVy = fopen(GetFilename("stokesVy"),"w");
    FILE * darcyVx  = fopen(GetFilename("darcyVx"),"w");
    FILE * darcyVy  = fopen(GetFilename("darcyVy"),"w");

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        fprintf(stokesVx, "%f ", ux[j*M_+i]);
        fprintf(stokesVy, "%f ", uy[j*M_+i]);
        fprintf(darcyVx, "%f ", vx[j*M_+i]);
        fprintf(darcyVy, "%f ", vy[j*M_+i]);
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
                        ux, uy, *br_, *basis_, parameter);

    CGNSPrepareParallel(&destDarcy_sol, &destDarcy_g, refArrayDarcyEssen_, mi,
                        vx, vy, *hdiv_, *basis_, parameter);

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

int Driver::PrintFlowEvent(){

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

    CGNSPrepareParallel(&sresult_->vel_stokes, &sresult_->g_stokes, 
                        refArrayStokesEssen_, mi, ux, uy, *br_, *basis_, parameter);   
    CGNSPrepareParallel(&sresult_->vel_darcy, &sresult_->g_darcy,
                        refArrayDarcyEssen_, mi, vx, vy, *hdiv_, *basis_, parameter);

//#ifdef CGNS_OUT
//    char stokesfile[] = "stokes.cgns";   
//    CgnsArrayOutput(dmMesh,&globalmesh,ux,uy,mi.MPIlocalCellStart[0],
//                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
//                    mi.MPIlocalCellSize[1],stokesfile);

//    char darcyfile[] = "darcy.cgns";    	
//    CgnsArrayOutput(dmMesh,&globalmesh,vx,vy,mi.MPIlocalCellStart[0],
//                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
//                    mi.MPIlocalCellSize[1],darcyfile);
//#endif

#ifndef CGNS_OUT
    quiverOutputEvent(ux,uy,vx,vy);
#endif

    return 1;
}

int Driver::PrintPorosity(){

    // Cell centered grid
    // Print porosity for the first time
    // No mesh files required
    // Compute porosity value with initial distribution of HD and CD
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

int Driver::PrintGrid(){

    FILE *gridPorox = fopen("gridCellX.dat", "w");
    FILE *gridPoroy = fopen("gridCellY.dat", "w");

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        vertex local {0.0,0.0};

        basis_->GetCorners(mi,{i,j});

        vertex global = GaussMapPointsFace(local, basis_->corners());

        if (withUnit){
            global[0] *= myPhase->pp->l0;
            global[1] *= myPhase->pp->l0;
        }

        fprintf(gridPorox,"%f ",global[0]);
        fprintf(gridPoroy,"%f ",global[1]);

    }
    fprintf(gridPorox, "\n");
    fprintf(gridPoroy, "\n");}

    fclose(gridPorox);
    fclose(gridPoroy);

    return 1;
}

int Driver::PrintPhaseEvent(){
    // Print porosity distribution when mesh file already exists
    // Use file passed by parameter
    // Simple print function for serial code, not parallel compatible

    //if (exists_file("girdCellX.dat") != 0 || exists_file("gridCellY.dat") != 0) {cout << "mesh file not exist!" << endl; return 0;}

    FILE *fp = fopen(GetFilename("porosity"), "w");
    FILE *ft = fopen(GetFilename("temperature"), "w");

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        vertex local {0.0,0.0};

        basis_->GetCorners(mi,{i,j});

        vertex global = GaussMapPointsFace(local, basis_->corners());
 
        double depth  = myPhase->pPtr->GetDepth(global[1], myPhase->pp->l0);  
        double lithoP = myPhase->pPtr->GetScaledLithoP(global[1]*(-1)*myPhase->pp->l0*0.6); 

        // Extract HD and CD from global solution vectors
        double HD, CD; // Get cell averaged values for approximation
        const int idx = FlatIndic(M_, {i,j});
        PetscCall(VecGetValues(globalHD, 1, &idx, &HD));
        PetscCall(VecGetValues(globalCD, 1, &idx, &CD));

        myPhase->pPtr->evalPhase(HD, CD, lithoP);

        //fprintf(fp, "%f ", myPhase->pPtr->phi.mlt);

        // Test ========================================================

        double temp = AssignPorosity(global, myPhase->pp);
        fprintf(fp, "%f ", temp);

        // =============================================================

        fprintf(ft, "%f ", myPhase->pPtr->TD);
    }
    fprintf(fp, "\n"); 
    fprintf(ft, "\n");}

    fclose(fp);
    fclose(ft);

    return 1;
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

int Driver::PrintPressureConstantOriginal(){

    // Print original pressure from dimensionless variables
    // p = rho_r g l_0 q + rho_f g z
   
    // Print piece wise constant pressure
    if (exists_file("gridCellX.dat") == 0 || exists_file("gridCellY.dat") == 0){
        // Checking for cell centered grid file
        // These two grid files are output along with printflow 
        cout << "No Grid File " << endl;
    } 

    FILE * darcypori  = fopen("darcypori.dat","w");
    FILE * stokespori = fopen("stokespori.dat","w");

    FILE * referencep = fopen("referencep.dat","w");

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

        vertex local {0.0,0.0};

        basis_->GetCorners(mi,{i,j});

        vertex globalver = GaussMapPointsFace(local, basis_->corners());

        // Retrieve original variables
        dp = myPhase->pp->rho_s * myPhase->pp->gy * myPhase->pp->l0*(-1*dp+globalver[1]);
        sp = myPhase->pp->rho_s * myPhase->pp->gy * myPhase->pp->l0*(-1*sp+globalver[1]);

        double rp = myPhase->pp->rho_s * myPhase->pp->gy * myPhase->pp->l0*(globalver[1]);

        fprintf(darcypori, "%f ", dp);
        fprintf(stokespori, "%f ", sp);
        fprintf(referencep, "%f ", rp);

    }fprintf(darcypori, "\n");
     fprintf(stokespori, "\n");
     fprintf(referencep, "\n");}

    fclose(stokespori);
    fclose(darcypori);
    fclose(referencep);

    return 1;
}

int Driver::PrintLithoPressure(){

    FILE * referencep = fopen("referencep.dat","w");

    int istart = mi.MPIlocalCellStart[0];
    int jstart = mi.MPIlocalCellStart[1];

    for (int j=jstart; j<jstart + mi.MPIlocalCellSize[1]; j++){
    for (int i=istart; i<istart + mi.MPIlocalCellSize[0]; i++){

        vertex local {0.0,0.0};

        basis_->GetCorners(mi,{i,j});

        vertex globalver = GaussMapPointsFace(local, basis_->corners());

        double rp = myPhase->pp->rho_s * myPhase->pp->gy * myPhase->pp->l0*(globalver[1]);

        fprintf(referencep, "%f ", rp);

    }fprintf(referencep, "\n");}

    fclose(referencep);

    return 1;
}

int Driver::PrintPressureEvent(){

    // Print pressure and pressure potentials according to events
    FILE * stokesq  = fopen(GetFilename("stokesq"),"w");
    FILE * darcyq = fopen(GetFilename("darcyq"),"w");
    FILE * stokesp  = fopen(GetFilename("stokesp"),"w");
    FILE * darcyp = fopen(GetFilename("darcyp"),"w");

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

        fprintf(stokesq, "%f ", sp);
        fprintf(darcyq, "%f ", dp);

        vertex local {0.0,0.0};

        basis_->GetCorners(mi,{i,j});

        vertex globalver = GaussMapPointsFace(local, basis_->corners());

        // Retrieve original variables
        dp = myPhase->pp->rho_s * myPhase->pp->gy * myPhase->pp->l0*(-1*dp+globalver[1]);
        sp = myPhase->pp->rho_s * myPhase->pp->gy * myPhase->pp->l0*(-1*sp+globalver[1]);

        double rp = myPhase->pp->rho_s * myPhase->pp->gy * myPhase->pp->l0*(globalver[1]);

        fprintf(stokesp, "%f ", sp);
        fprintf(darcyp, "%f ", dp);

    }fprintf(stokesq, "\n");
     fprintf(darcyq, "\n");
     fprintf(stokesp, "\n");
     fprintf(darcyp, "\n");}

    fclose(stokesq);
    fclose(darcyq);
    fclose(stokesp);
    fclose(darcyp);

    return 1;
}

int Driver::PrintBoundaryDOFs(){

    // Prevent mesh from being too large
	 // Do a upto 3*3 mesh 
	 // Print dof information to screen with restriction to mesh size
    assert(M_ < 3);
    assert(N_ < 3);

    std::cout << "Total number of essential dof for Stokes : " << bndryDOFStokes_ << std::endl;

    for (int i=0; i<br_->getDOF(); i++){
        // Print out all dirichlet dofs
        std::cout << "Global dof : " << i << " Actual dof pos: "  << 
                      refArrayStokesEssen_[i] << std::endl;
    }

    std::cout << "Print Natural boundary condition in details."  << std::endl;


    for (const auto & [key, value] : refArrayStokesNatur_){
        std::cout << key << " : " << value << std::endl;
    }

    std::cout << std::endl;

/*
    for (int i=0; i<hdiv_->getDOF(); i++){
        std::cout << "Dirichlet dofs for Darcy: " << std::endl;
        std::cout << refArrayDarcyEssen_[i] << std::endl;
    }
*/

    return 1;
}

int Driver::PrintStokesBoundaryDOFs(){

    FILE * stokesBndryDOF = fopen("stokesBndryDOF.txt","w");

    fprintf(stokesBndryDOF, "Boundary dof details for stokes: \nEssential boundary dofs : %d \nNatural boundary dofs : %d \n\n", bndryDOFStokes_, bndryDOFStokesNatur_);

    fprintf(stokesBndryDOF, "Details of essential boundary dofs: \n"); 

    //cout << br_->getDOF() << endl;

    // Print all essential boundary dofs to the file
    for (int i=0; i<br_->getDOF(); i++){
        std::vector<int> localinfo = br_->GlobalToLocalMapBndry(mi, i);
        if (localinfo.at(0) > -1){
            //Interior dof will be marked -1
            bndryType type = bndryTypeMarker(mi, Bend(mi, localinfo.at(0)),localinfo.at(1), {0}); 

            fprintf(stokesBndryDOF, "Global dof index : %d, global cell index : %d, local dof index : %d ", i, localinfo.at(0), localinfo.at(1));
            if (type == dirichlet){
                fprintf(stokesBndryDOF, " Dirichlet\n");
            } else {
                fprintf(stokesBndryDOF, " Neumann \n");
            }
        }
    }

    fclose(stokesBndryDOF);

    return 1;
}
