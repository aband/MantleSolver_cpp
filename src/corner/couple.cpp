// Another coupling method
#include "couple.h"

int couple::CreatePhase(){

    myPhase = Phase();
    myPhase.pp = (PhysProperty *)malloc(sizeof(PhysProperty));

    AssignPhyProperties(myPhase.pp);

    myPhase.pPtr = new EUTECTIC::phase();

    return 1;
}

int couple::ShowPhase(){

    // Showing phase attributes
    cout << " ========================================================= " << endl;
    cout << "Phase attributes defined in eutectic phase class ...       " << endl;

    myPhase.pPtr->printInfo();

    cout << " ========================================================= " << endl;
    cout << "Phase attributes defined in AssignPhyProperties function .." << endl;
    cout << "Compaction length        : " << myPhase.pp->l0 << " m" << endl;
    cout << "Upwelling solid velocity : " << myPhase.pp->V0 <<" m/s, " << 
            myPhase.pp->V0*365*24*3600*100 << " cm/yrs "<< endl;
    cout << "Characteristic velocity  : " << -1 *myPhase.pp->u0 << " m/s" << endl;
    cout << "characteristic time step : " << abs(myPhase.pp->l0 / myPhase.pp->u0) << " s , " 
                                          << abs(myPhase.pp->l0/myPhase.pp->u0 /365/24/3600) << " yrs"<< endl;
    cout << "Characteristic permeability: " << 1.0/myPhase.pp->invk0 << " m^2" << endl;
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
        double physscale = myPhase.pp->L0/myPhase.pp->l0;
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

    N_ = N;
    M_ = M;

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
    // Calculate total dofs    
    // Edge dofs are always vertical edges counted first
 
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    const valarray<double>& gwf = GaussWeightsFace; 
    const vector<vertex>& gpf = GaussPointsFace;

    // ==================================================

    int tolvertgauss = (M+1)*N*gpe.size();
    int tolhorigauss = M*(N+1)*gpe.size();

    int toledgegauss = tolvertgauss + tolhorigauss;

    edgegauss.resize(toledgegauss);

    int tolcellgauss = M*N*gwf.size();

    cellgauss.resize(tolcellgauss);

    cellcenter.resize(M*N);

    vertexSet vertedge;
    vertexSet horiedge;

    int indexvert = 0;
    int indexhori = 0;

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        indice gcell {i,j};
        vertexSet corners = extractCorners(mi, gcell);

        vertedge = {corners.at(0), corners.at(3)};
        horiedge = {corners.at(0), corners.at(1)};

        indexvert = (j*(M+1) + i)*gpe.size();
        indexhori = tolvertgauss + (j*M + i)*gpe.size();

        for (int g=0; g<gpe.size(); g++){
            edgegauss.at(indexvert + g) = GaussMapPointsEdge({gpe[g]},vertedge);
            edgegauss.at(indexhori + g) = GaussMapPointsEdge({gpe[g]},horiedge);
        }   

        // cell centered and cell gauss points
        cellcenter.at(j*M+i) = (corners.at(0) + corners.at(1) + 
                                corners.at(2) + corners.at(3))/4.0;

        for (int g=0; g<gwf.size(); g++){
            cellgauss.at(gwf.size()*(j*M+i) + g) = 
                         GaussMapPointsFace(gpf[g],corners); 
        }

    }}

    // right vertical edge
    for (int j=0; j<N; j++){

        indice gcell {M-1, j};
        vertexSet corners = extractCorners(mi, gcell);

        vertedge = {corners.at(1), corners.at(2)};

        indexvert = (j*(M+1) + M)*gpe.size();
        for (int g=0; g<gpe.size(); g++){
            edgegauss.at(indexvert + g) = GaussMapPointsEdge({gpe[g]}, vertedge);
        }

    }

    // Top horizontal edge
    for (int i=0; i<M; i++){

        indice gcell {i, N-1};
        vertexSet corners = extractCorners(mi, gcell);

        horiedge = {corners.at(3), corners.at(2)};

        indexhori = tolvertgauss + (N*M + i)*gpe.size();

        for (int g=0; g<gpe.size(); g++){
            edgegauss.at(indexhori + g) = GaussMapPointsEdge({gpe[g]}, horiedge);
        }

    }

    return 1;
}

//int couple::PrepareFlow(){

    //ds = DarcyStokes();

    //ds.init(mi, myPhase->pp, {0.0});

    //ds.Assemble(mi, edgeporo, cellporo, average_poro, {0.0});

    //return 1;
//}
