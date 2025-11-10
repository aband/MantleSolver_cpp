// Another coupling method
#include "couple.h"

int couple::CreatePhase(){

    myPhase = new Phase();
    myPhase->pp = (PhysProperty *)malloc(sizeof(PhysProperty));

    AssignPhyProperties(myPhase->pp);

    myPhase->pPtr = new EUTECTIC::phase();

    return 1;
}

int couple::ShowPhase(){

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

int couple::CreateMesh(const int& M, const int& N,
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


    // Compute and store all the gauss points



    return 1;
}
