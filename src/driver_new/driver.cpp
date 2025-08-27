#include "driver.h"

int Driver::CreatePhase(){

    myPhase = new Phase();

    myPhase->pp = (PhysProperty *)malloc(sizeof(PhysProperty));

    AssignPhyProperties(myPhase->pp);

    myPhase->pPtr = new EUTECTIC::phase();

    return 0;
}

int Driver::ShowPhase(){

    // Showing phase attributes
    cout << " ========================================================= " << endl;
    cout << "Phase attributes defined in eutectic phase class ...       " << endl;

    myPhase->pPtr->printInfo();

    cout << " ========================================================= " << endl;
    cout << "Phase attributes defined in AssignPhyProperties function .." << endl;
    cout << "Compaction length        : " << myPhase->pp->l0 << " m" << endl;
    cout << "Upwelling solid velocity : " << myPhase->pp->V0 <<" m/s, " << 
            myPhase->pp->V0*365*24*3600*100 << " cm/yrs "<< endl;
    cout << "Characteristic velocity  : " << -1 *myPhase->pp->u0 << " m/s" << endl;
    cout << "characteristic time step : " << abs(myPhase->pp->l0 / myPhase->pp->u0) << " s , " 
                                          << abs(myPhase->pp->l0/myPhase->pp->u0 /365/24/3600) << " yrs"<< endl;
    cout << "Characteristic permeability: " << 1.0/myPhase->pp->invk0 << " m^2" << endl;
    cout << "Scaled characteristic permeability: "     << endl;
    cout << " ========================================================= " << endl;

    return 1;
}

int Driver::CreateMesh(const int& M, const int& N,
                       double L, double H, 
                       double xstart, double ystart,
                       const int& stencilWidthMesh, 
                       const int& stencilWidthU,
                       const bool& physicsScale,
                       const int& meshType){

    if (physicsScale){
        double physscale = myPhase->pp->L0/myPhase->pp->l0;
        L = L*physscale;
        H = H*physscale;
        xstart = xstart*physscale, 
        ystart = ystart*physscale;
    }

    // Create dmMesh
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidthMesh, NULL, NULL, 
    &dmMesh));
    PetscCall(DMSetFromOptions(dmMesh));              
    PetscCall(DMSetUp(dmMesh));

    // Create dmU
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 1, 
    stencilWidthU, NULL, NULL, &dmu));
    PetscCall(DMSetFromOptions(dmu));              
    PetscCall(DMSetUp(dmu));     

    // Create MeshParam object (historical object one time use only)
    MeshParam mp;
    mp.xstart = xstart;
    mp.ystart = ystart;
    mp.L = L;
    mp.H = H;

    mi.L = L;
    mi.H = H;
    L_ = L;
    H_ = H;

    // Create global vector containing mesh
    PetscCall(DMCreateGlobalVector(dmMesh, &globalmesh));
    switch(meshType){
        case 0: CreateFullMesh(dmMesh, &globalmesh, &mp); break;
        case 1: LogicRectMesh(dmMesh, &globalmesh, &mp);  break;
        case 2: RefineMesh(dmMesh, &globalmesh, &mp);
        //case 2: TestControlMeshSecond(dmCell,L,H); break;
        //case 3: TestControlMeshThird(dmCell,L,H);  break;
    }

    ReadMeshPortion(dmMesh, &globalmesh, mi.lmesh);

    AssignValuesMeshInfo(mi, dmMesh, dmu);

    return 1;
}

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

    CreateRefMap(*br_  , mi, refArrayStokesEssen_, refArrayStokesNatur_, &bndryDOFStokes_, &bndryDOFStokesNatur_, parameter);
    CreateRefMap(*hdiv_, mi, refArrayDarcyEssen_ , refArrayDarcyNatur_ , &bndryDOFDarcy_ , &bndryDOFDarcyNatur_ , parameter);

    Result_ = (ReducedSys *)malloc(sizeof(ReducedSys));

    sresult_ = (ScatterResult *)malloc(sizeof(ScatterResult));

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

int Driver::exactandreconstructTest(){

    // Eat and spit test
    printexactsol(mi, 0, InitCD, 1, true, {0.0});

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    vector<double> sigma_lg;
    sigma_lg.resize(stenlg.size());

    vector<double> sigma_sm;
    sigma_sm.resize(stensm.size());

    Vec localvecCD; 
    double ** locvalsCD;

    // Distribute local part to local vectors.
    PetscCall(DMGetLocalVector(dmu, &localvecCD)); 

    PetscCall(DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localvecCD));
    PetscCall(DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localvecCD));

    PetscCall(DMDAVecGetArray(dmu, localvecCD, &locvalsCD));

    for (int s=0; s<stenlg.size(); s++){
        sigma_lg.at(s) = stenlg.at(s).sigma(locvalsCD);
    }

    for (int s=0; s<stensm.size(); s++){
        sigma_sm.at(s) = stensm.at(s).sigma(locvalsCD);
    }

    // Setup nonlinear weights
    for (int s=0; s<my_recon_CD.size(); s++){
        my_recon_CD.at(s).extractsigma(sigma_lg, sigma_sm);
        my_recon_CD.at(s).setWgts(L_*H_/(double)M/(double)N);
    }

    printreconsol2(my_recon_CD, M, N, 1, stenlg, stensm, mi, locvalsCD);

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
				int s = j*M+i;
        cout << my_recon_CD.at(s).efforder() << "  ";
    }cout << endl;}

    DMDAVecRestoreArray(dmu,localvecCD,&locvalsCD);
    DMRestoreLocalVector(dmu, &localvecCD); 

    cout << endl << endl;

    // ==========================================================
    printexactsol(mi, 0, InitHD, 2, true, {0.0});

    sigma_lg.clear();
    sigma_lg.resize(stenlg.size());

    sigma_sm.clear();
    sigma_sm.resize(stensm.size());

    Vec localvecHD; 
    double ** locvalsHD;

    // Distribute local part to local vectors.
    PetscCall(DMGetLocalVector(dmu, &localvecHD)); 

    PetscCall(DMGlobalToLocalBegin(dmu, globalHD, INSERT_VALUES, localvecHD));
    PetscCall(DMGlobalToLocalEnd(dmu, globalHD, INSERT_VALUES, localvecHD));

    PetscCall(DMDAVecGetArray(dmu, localvecHD, &locvalsHD));

    for (int s=0; s<stenlg.size(); s++){
        sigma_lg.at(s) = stenlg.at(s).sigma(locvalsHD);
    }

    for (int s=0; s<stensm.size(); s++){
        sigma_sm.at(s) = stensm.at(s).sigma(locvalsHD);
    }

    // Setup nonlinear weights
    for (int s=0; s<my_recon_HD.size(); s++){
        my_recon_HD.at(s).extractsigma(sigma_lg, sigma_sm);
        my_recon_HD.at(s).setWgts(abs(L_*H_)/(double)M/(double)N);
    }

    printreconsol2(my_recon_HD, M, N, 2, stenlg, stensm, mi, locvalsHD);

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
				int s = j*M+i;
        cout <<"At cell " << i << ", " << j  << "  " << my_recon_HD.at(s).efforder() << "   ";
		  cout << my_recon_HD.at(s).printinfo() << endl;
    }cout << endl;}

    DMDAVecRestoreArray(dmu,localvecHD,&locvalsHD);
    DMRestoreLocalVector(dmu, &localvecHD); 

    return 1;
}
