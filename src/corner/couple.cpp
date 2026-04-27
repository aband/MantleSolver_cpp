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

    phasequantityedge.resize(toledgegauss);

    int tolcellgauss = M*N*gwf.size();

    cellgauss.resize(tolcellgauss);

    cellcenter.resize(M*N);

    phasequantitycell.resize(M*N);

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

int couple::PrepareTransport(TransportVariable& H,
                             TransportVariable& C,
                             double (*initHD)(const valarray<double>& point, 
                                              const vector<double>& param),
                             double (*initCD)(const valarray<double>& point, 
                                              const vector<double>& param)){

    // Create global solution vectors
    PetscCall(DMCreateGlobalVector(dmu, &H.sol)); 
    PetscCall(DMCreateGlobalVector(dmu, &C.sol)); 

    // Assign cell averaged values as initial condition
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &H.sol, {H_,0.0}, initHD);
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &C.sol, {H_,0.0}, initCD);

    return 1;
}

int couple::ReadVectorTransport(Vec * H, Vec * C, 
					                 const char * fileH,
										  const char * fileC,
										  int mark){

    Vec tH = *H;
    Vec tC = *C;

//    FILE * fC = fopen(GetFilename(fileC, mark), "r");
//    FILE * fH = fopen(GetFilename(fileH, mark), "r");
    FILE * fC = fopen("cellC1.dat", "r");
    FILE * fH = fopen("cellH1.dat", "r");

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        double valH = 0.0; 
		  double valC = 0.0;

        indice index {i,j};
		  int nelem = FlatIndic(mi, index);

        fscanf(fC, "%lf", &valC);
        fscanf(fH, "%lf", &valH);

        PetscCall(VecSetValue(tH, nelem, valH, INSERT_VALUES));
        PetscCall(VecSetValue(tC, nelem, valC, INSERT_VALUES));

    }}

    PetscCall(VecAssemblyBegin(tH));
    PetscCall(VecAssemblyEnd(tH));

    PetscCall(VecAssemblyBegin(tC));
    PetscCall(VecAssemblyEnd(tC));

	 fclose(fC);
	 fclose(fH);

    return 1;
}

/*
int couple::PhaseSplit(TransportVariable& H, 
                       TransportVariable& C){

    double thisH = 0.0;
    double thisC = 0.0;

    // Pressure
    double pressure = 0.0;

    int celldof = 0;
    int edgedof = 0;
    int cellgaussdof = 0; 

    // Clear all the phase variables 
    edgeporo.clear();
    edgeporo.resize(M_*N_*12);

    cellporo.clear();
    cellporo.resize(M_*N_*9);

    average_poro.clear(); 
    average_poro.resize(M_*N_);

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){
       
        celldof = j*M_+i;

        // Edge gauss points
        for (int e=0; e<3*4; e++){
            edgedof = celldof*12 + e;
            vertex edgep = edgegauss.at(edgedof);

            pressure = myPhase.pPtr->GetStaticP(-1*edgep[1], myPhase.pPtr->l0);

            thisH = H.my_recon.at(celldof)->elem_val.at(e);
            thisC = C.my_recon.at(celldof)->elem_val.at(e);

            myPhase.pPtr->evalPhase(thisH, thisC, pressure);

            // Store evaluated values
            edgeporo.at(edgedof) = myPhase.pPtr->pc.phil;
        }

        // Cell gauss points
        for (int g=0; g<9; g++){
            cellgaussdof = celldof*9 + g;
            vertex cellp = cellgauss.at(cellgaussdof);

            pressure = myPhase.pPtr->GetStaticP(-1*cellp[1], myPhase.pPtr->l0);

            thisH = H.cellgauss.at(cellgaussdof);
            thisC = C.cellgauss.at(cellgaussdof);

            myPhase.pPtr->evalPhase(thisH, thisC, pressure);
            cellporo.at(cellgaussdof) = myPhase.pPtr->pc.phil; 
        }

        // Center point
        vertex centerp = cellcenter.at(celldof);

        pressure = myPhase.pPtr->GetStaticP(-1*centerp[1], myPhase.pPtr->l0);

        thisH = H.cellcenter.at(celldof);
        thisC = C.cellcenter.at(celldof);
    }}

    return 1;
}
*/
int couple::ExtractThisEdgeGauss(vector<vertex>& thisedgegauss, int gsize,
					                  int i, int j, int xSize, int ySize, int offset){

   int dof = j*xSize + i + offset;

   for (int g=0; g<gsize; g++){
       thisedgegauss.at(g) = edgegauss.at(dof*gsize + g);
	}

    return 1;
}

static int combineEdgePhase(EUTECTIC::PhaseComp& pcneg, 
                            EUTECTIC::PhaseComp& pcpos, 
				                EUTECTIC::PhaseComp& pcedge){

    // geometric average of phi
	 if (pcneg.phil * pcpos.phil == 0.0){
        pcedge.phil == 0.0;
	 }	else {
        pcedge.phil = harmonic_mean(pcneg.phil, pcneg.phil);
	 }

    pcedge.cl = (pcneg.cl + pcpos.cl)/2.0;
    pcedge.cs = (pcneg.cs + pcpos.cs)/2.0;

    pcedge.TDp = (pcneg.TDp + pcpos.TDp)/2.0;

    return 1;
}

int couple::EvaluateThisEdgePhase(int gsize, int i, int j, int xSize, int ySize, int offset, 
								          const vector<double>& Hneg, const vector<double>& Cneg, 
								          const vector<double>& Hpos, const vector<double>& Cpos){

    int dof = j*xSize + i + offset;

	 EUTECTIC::PhaseComp pcpos, pcneg;

    for (int g=0; g<gsize; g++){
		  vertex gaussp = edgegauss.at(dof*gsize + g);
        double lithoP = myPhase.pPtr->GetStaticP(-1*gaussp[1],
                                                  myPhase.pPtr->l0);

        myPhase.pPtr->evalPhase(Hneg.at(g), Cneg.at(g), lithoP);
        myPhase.pPtr->copyPhaseComp(pcneg);

        myPhase.pPtr->evalPhase(Hpos.at(g), Cpos.at(g), lithoP);
        myPhase.pPtr->copyPhaseComp(pcpos);

        combineEdgePhase(pcneg, pcpos, phasequantityedge.at(dof*gsize+g));
	 }

    return 1;
}

int couple::computePorosityEdgePhase(const MeshInfo& mi,
					                     int xSize, int ySize, int offset,
                                    int xMaxCell, int yMaxCell,
                                    TransportVariable& H,
												TransportVariable& C){

  	 const valarray<double>& gpe = GaussPointsEdge;

    bool onbndry = false;
    int bndrytype = 0;

    int edgepos = 0;
    int edgeneg = 0;

    indice cellpos {0,0};
    indice cellneg {0,0};

    // Reconstruction values evaluated from neg cell
    vector<double> Hneg;
	 vector<double> Cneg;
    // Reconstruction values evaluated from pos cell
    vector<double> Hpos;
	 vector<double> Cpos;

    vector<vertex> thisedgegauss;
    thisedgegauss.resize(gpe.size());

    for (int j=0; j<ySize; j++){
    for (int i=0; i<xSize; i++){
        H.getNeighbors(xMaxCell, yMaxCell, xSize, ySize, i, j, 
			              edgepos, edgeneg, cellpos, cellneg, onbndry);

        int dof = j*xSize + i;
        // Edge dof
        dof += offset;

        // Extract velocity current edge
        H.ExtractThisEdge(mi, cellneg, edgeneg, 
                              cellpos, edgepos, gpe.size(),
                              Hneg, Hpos);

        C.ExtractThisEdge(mi, cellneg, edgeneg, 
                              cellpos, edgepos, gpe.size(),
                              Cneg, Cpos);

        // Extract edge gauss points
        //ExtractThisEdgeGauss(thisedgegauss, gpe.size(), i, j, xSize, ySize, offset); 
        EvaluateThisEdgePhase(gpe.size(), i, j, xSize, ySize, offset, Hneg, Cneg, Hpos, Cpos);

	 }}

    return 1;
}

// Compute cell averaged value can be replaced
int couple::computePorosityCellPhase(TransportVariable& H, 
												 TransportVariable& C){

    double cellH, cellC;

    average_poro.clear();
    average_poro.resize(M_*N_);

    cellphase.clear();
    cellphase.resize(M_*N_);

    meltp.clear();
    meltp.resize(M_*N_);

    currentTemp.clear();
    currentTemp.resize(M_*N_);

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        int dof = j*M_ + i;

        VecGetValues(H.sol, 1, &dof, &cellH);
        VecGetValues(C.sol, 1, &dof, &cellC);

        vertex center = cellcenter.at(dof);

        double lithoP = myPhase.pPtr->GetStaticP(-1*center[1],
                                                  myPhase.pPtr->l0);

        myPhase.pPtr->evalPhase(cellH, cellC, lithoP);

        average_poro.at(dof) = myPhase.pPtr->pc.phil;

		  cellphase.at(dof) = myPhase.pPtr->pc.region;

		  currentTemp.at(dof) = myPhase.pPtr->pc.TDp;

        meltp.at(dof) = myPhase.pPtr->pc.Tep;

	 }}

    return 1;
}

int couple::computePorosity_phase(TransportVariable& H,
                                  TransportVariable& C){

    const valarray<double>& gpe = GaussPointsEdge;

    edgeporo.clear();   
    edgeporo.resize(edgegauss.size());

    // Edge phase quantities, with harmonic average
    computePorosityEdgePhase(mi, M_+1, N_,         0, M_, N_, H, C);
    computePorosityEdgePhase(mi, M_, N_+1, (M_+1)*N_, M_, N_, H, C);

    for (int i=0; i<(int)edgeporo.size(); i++){

        edgeporo.at(i) = phasequantityedge.at(i).phil;

	 }

    cellporo.clear();
    cellporo.resize(M_*N_*9);

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){
        int dof = j*M_ + i; 
        for (int g=0; g<9; g++){

            int gdof = dof*9 + g;
            vertex point = cellgauss.at(gdof); 
            double lithop = myPhase.pPtr->GetStaticP(-1*point[1], myPhase.pPtr->l0);

            double thisH = H.cellgauss.at(gdof);
            double thisC = C.cellgauss.at(gdof);

            myPhase.pPtr->evalPhase(thisH, thisC, lithop);
            cellporo.at(gdof) = myPhase.pPtr->pc.phil;

        }
    }}

    // Cell centerred phase quantities
    computePorosityCellPhase(H, C); 

    return 1;
}

int couple::calculatePhaseVel(const vector<vertex>& stokesvel,
					               const vector<vertex>& relativevel){

	 // Resize velocity vectors
    effvel.resize(stokesvel.size());
    phasevel.resize(stokesvel.size());
    solidvel.resize(stokesvel.size());

    for (int g=0; g<(int)stokesvel.size(); g++){

        double cl   = phasequantityedge.at(g).cl;
        double cs   = phasequantityedge.at(g).cs;
        double phil = phasequantityedge.at(g).phil;

        effvel.at(g) = cl*phil* (stokesvel.at(g) + relativevel.at(g)) + cs*(1-phil)*stokesvel.at(g);
        effvel.at(g) /= cl*phil + cs*(1-phil);

        phasevel.at(g) = phil*relativevel.at(g) + stokesvel.at(g);

        solidvel.at(g) = (1-phil) * stokesvel.at(g);
	 }

    return 1;
}
