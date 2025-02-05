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

int Driver::PrepareTransport(double (*funcHD)(const valarray<double>& point, 
                                              const vector<double>& param),
                             double (*funcCD)(const valarray<double>& point, 
                                              const vector<double>& param)){

    PetscCall(DMCreateGlobalVector(dmu, &globalCD));
    PetscCall(DMCreateGlobalVector(dmu, &globalHD));

    // Assign cell averaged values as initial condition
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalCD, {H_,0.0}, funcCD);
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalHD, {myPhase->pp->l0*H_,0.0}, funcHD);

    // Initialization of multi level weno and corresponding usage
    ml = multilevel(); 

    ml.addLevel("(3,3)", {3,3}, mi);
    ml.addLevel("(2,2)", {2,2}, mi);

    advection = mluse();

    // Test for nonlinear weighting
    unordered_map<std::string, vector<indice>> method;
    method.insert(std::make_pair<std::string, vector<indice>>("(3,3)", { {-1,-1} }));
    method.insert(std::make_pair<std::string, vector<indice>>("(2,2)", { {-1,-1}, {0,-1}, {0,0}, {-1,0} }));

    // Area scale
    h0 = sqrt((L_*H_)/
         (double)(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]));

    advection.setmethod("all", method);
    advection.setbias("all");

    // Compute bottom fixed value
    HDbottom = funcHD({0.0,-1*H_},{myPhase->pp->l0*H_,0.0});
    CDbottom = 0.1;

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

int Driver::SolveFlow(int maxIter, double tolUzawa, const Tensor<weights>& allwgtsHD, double ** lHD, 
                                                    const Tensor<weights>& allwgtsCD, double ** lCD){

    ParallelMatrixAssemble(allwgtsHD, lHD, allwgtsCD, lCD);

    int nelem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    CreateLinearSys(reducedStokes_, nelem);
    CreateLinearSys(reducedDarcy_, nelem);

    CreateCoupledSystem(reducedStokes_, reducedDarcy_, Result_, &K);

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
