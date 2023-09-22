#include "driver.h"

Driver::~Driver(){

    // Clear used objects
    DMDAVecRestoreArray(dmu_,localu_,&mi_.localVals);
    DMRestoreLocalVector(dmu_, &localu_); 

    VecDestroy(&globalu_); 
    VecDestroy(&fullmesh_);
    DMDestroy(&dmu_);
    DMDestroy(&dmMesh_);
}

int Driver::Prepare(){

    PetscErrorCode ierr;

    // Default global cell size
    globalM_ = 3;
    globalN_ = 3;

    PetscCall(PetscOptionsGetInt(NULL,NULL,"-M",&globalM_,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-N",&globalN_,NULL));

    // Create data management object for solution.
    int stencilWidthMesh_ = 5;

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    globalM_, globalN_, PETSC_DECIDE, PETSC_DECIDE, 2, 
    stencilWidthMesh_, NULL, NULL, &dmMesh_));       
    PetscCall(DMSetFromOptions(dmMesh_));              
    PetscCall(DMSetUp(dmMesh_));
    PetscCall(DMCreateGlobalVector(dmMesh_, &fullmesh_));

    // Define physical domain
    // The default size is from -1 to 1
    L_ = 2.0, H_ = 2.0;
    double xstart_ = -1.0, ystart_ = -1.0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L_,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H_,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart_, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart_, NULL));

    int singleStencilTest_ = 0;
    double scale_ = 1;
    PetscCall(PetscOptionsGetInt(NULL,NULL, "-single", &singleStencilTest_, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL, "-scale", &scale_, NULL));

    if (singleStencilTest_){
        L_ = L_/scale_;
        H_ = H_/scale_;
        xstart_ = -L_/2.0;
        ystart_ = -H_/2.0;
    }

    // Define different types of mesh.
    MeshParam mp;
    mp.xstart = xstart_;
    mp.ystart = ystart_;
    mp.L = L_;
    mp.H = H_;

    int meshtype_ = 0; 

    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshtype_,NULL));
    switch(meshtype_){
        case 0: CreateFullMesh(dmMesh_, &fullmesh_, &mp); break;
        case 1: LogicRectMesh(dmMesh_, &fullmesh_, &mp);  break;
        //case 2: TestControlMeshSecond(dmCell,L,H); break;
        //case 3: TestControlMeshThird(dmCell,L,H);  break;
    }

    int printmesh=0;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-printmesh",&printmesh,NULL));
    if(printmesh){ 
        VecView(fullmesh_, PETSC_VIEWER_STDOUT_WORLD);
        PrintFullMesh(dmMesh_, &fullmesh_);
    }

    // Create data management for solution
    int stencilWidthU_ = 3;

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    globalM_, globalN_, PETSC_DECIDE, PETSC_DECIDE, 1, 
    stencilWidthU_, NULL, NULL, &dmu_));
    PetscCall(DMSetFromOptions(dmu_));              
    PetscCall(DMSetUp(dmu_));                       

    // Create global vector
    PetscCall(DMCreateGlobalVector(dmu_,&globalu_));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Initialization finished ...\n"));

    return 0;
}

void Driver::CellAveragedInit(double (*func)(const valarray<double>& point,
                                             const vector<double>& param)){

    SimpleInitialValue(dmMesh_, dmu_, &fullmesh_, &globalu_, {-L_/(double)(2*globalM_)},func);
}

int Driver::CreateMeshInfo(){

    ReadMeshPortion(dmMesh_, &fullmesh_, mi_.lmesh);

    // Create local vector
    PetscCall(DMGetLocalVector(dmu_, &localu_));

    PetscCall(DMGlobalToLocalBegin(dmu_, globalu_, INSERT_VALUES, localu_));
    PetscCall(DMGlobalToLocalEnd(dmu_, globalu_, INSERT_VALUES, localu_));

    DMDAVecGetArray(dmu_, localu_, &mi_.localVals);
  
    AssignValuesMeshInfo(mi_, dmMesh_, dmu_);

    return 0;
}
