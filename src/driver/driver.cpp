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

    cout << "Compaction Length: " << myPhase->pp->l0 << " m" << endl;
    cout << "Upwelling solid velocity: " << myPhase->pp->V0 << " m/s" << endl;
    cout << "Characteristic Velocity: " << myPhase->pp->u0 << " m/s" << endl;
    cout << "Time step: " << endl;
    cout << "Characteristic permeability: " << endl;
    cout << "Characteristic Enthalpy: "     << endl;

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
    M, N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidthMesh, NULL, NULL, &dmMesh));
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
    mp_.xstart = xstart;
    mp_.ystart = ystart;
    mp_.L = L;
    mp_.H = H;

    // Create global vector containing mesh
    PetscCall(DMCreateGlobalVector(dmMesh, &globalmesh));

    switch(meshType){
        case 0: CreateFullMesh(dmMesh, &globalmesh, &mp_); break;
        case 1: LogicRectMesh(dmMesh, &globalmesh, &mp_);  break;
        case 2: RefineMesh(dmMesh, &globalmesh, &mp_);
        //case 2: TestControlMeshSecond(dmCell,L,H); break;
        //case 3: TestControlMeshThird(dmCell,L,H);  break;
    }

    ReadMeshPortion(dmMesh, &globalmesh, mi.lmesh);

    M_ = M;
    N_ = N;
  
    return 0;
}

int Driver::PrintMesh(){

    VecView(globalmesh, PETSC_VIEWER_STDOUT_WORLD);
    PrintFullMesh(dmMesh, &globalmesh);

    return 0;
}

char * Driver::GetFilename(const char * fieldname){

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];
    std::sprintf(n_char,"%d",eventCount);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    return filename;
}

int Driver::InitTransport(double (*funcHD)(const valarray<double>& point, const vector<double>& param),
                          double (*funcCD)(const valarray<double>& point, const vector<double>& param)){

    eventCount = 0;

    // Create global vectors
    PetscCall(DMCreateGlobalVector(dmu,&globalHD));
    PetscCall(DMCreateGlobalVector(dmu,&globalCD));

    // Assign Initial values in the form of cell-averaged value
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalCD, {0.0,0.0}, funcCD); 
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalHD, {myPhase->pp->l0,0.0}, funcHD); 

    // Distribute global to local vectors
    DMGetLocalVector(dmu, &localHD);

    DMGlobalToLocalBegin(dmu, globalHD, INSERT_VALUES, localHD);
    DMGlobalToLocalEnd(dmu, globalHD, INSERT_VALUES, localHD);

    DMGetLocalVector(dmu, &localCD);

    DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localCD);
    DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localCD);

    DMDAVecGetArray(dmu, localCD, &mi.localCD);
    DMDAVecGetArray(dmu, localHD, &mi.localHD);

    AssignValuesMeshInfo(mi, dmMesh, dmu);

    // Preparation for MLWENO
    mlpPtr_ = new MLWENO::MLWENOPrepare();

    // Initialize MLWENO objects
    mluseAdv_ = new MLWENO::MLWENOUse();

    mluseDif_ = new MLWENO::MLWENOUse();

    return 0;
}

int Driver::clean(){

    // Restore local vectors
    DMDAVecRestoreArray(dmu,localCD,&mi.localCD);
    DMRestoreLocalVector(dmu, &localCD); 
    DMDAVecRestoreArray(dmu,localHD,&mi.localHD);
    DMRestoreLocalVector(dmu, &localHD); 

    PetscCall(VecDestroy(&globalHD));
    PetscCall(VecDestroy(&globalCD));
 
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    return 0;
}
