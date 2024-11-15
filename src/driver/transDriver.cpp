#include "driver.h"

int Driver::AddLevels(const int& stencilSize){

    mlpPtr_->AddLevel(mi, stencilSize, stencilSize);

    return 1;
}

int Driver::AddLevels(const int& stencilSizeX,
                      const int& stencilSizeY){

    mlpPtr_->AddLevel(mi,stencilSizeX, stencilSizeY);

    return 1;
}

int Driver::AddLevels(const vector<int>& stencilSizes){


    for (const auto& stencilSize: stencilSizes){
        AddLevels(stencilSize);
    }

    return 1;
}

int Driver::AddLevels(const vector<pair<int, int>>& stencilSizes){


    for (const auto& stencilSize: stencilSizes){
        AddLevels(stencilSize.first, stencilSize.second);
    }

    return 1;
}

int Driver::PrepareDefaultTransport(){

    // Prepare levels 
    AddLevels(1);
    AddLevels(2);
    AddLevels(3);

    AddLevels(4,3);
    AddLevels(3,4);

    // Create set containing position information and field information
    posSet_.insert("interior");
    posSet_.insert("corner");
    posSet_.insert("edge");

    locFuncSet_["edge"] = edge; 
    locFuncSet_["corner"] = corner; 
    locFuncSet_["interior"] = interior; 

    fieldSet_.insert("HD");
    fieldSet_.insert("CD");

    // Create Smoothness Indicator for the first time
    mlpPtr_->UpdateSmoothnessIndic(mi, mi.localCD, "CD");
    mlpPtr_->UpdateSmoothnessIndic(mi, mi.localHD, "HD");

    // Three different treatment on interior, edge and corner cells
    mluseAdv_->AddMLWENOLevel("interior", {"(3,3)","(2,2)"}, mlpPtr_);

    mluseAdv_->AssignWENOStencils("interior","(2,2)",{{-1,0},{-1,-1},{0,0},{0,-1}});
    mluseAdv_->AssignWENOStencils("interior","(3,3)",{{-1,-1}});
    mluseAdv_->AssignLinearWgts("interior","(2,2)",{1,1,1,1});
    mluseAdv_->AssignLinearWgts("interior","(3,3)",{5});
    mluseAdv_->UpdateNonLinearWgts(mi, "interior", interior, "HD"); 
    mluseAdv_->UpdateNonLinearWgts(mi, "interior", interior, "CD"); 

    mluseAdv_->AddMLWENOLevel("edge", {"(3,3)", "(2,2)"}, mlpPtr_);
    mluseAdv_->AssignWENOStencils("edge","(2,2)",{{-1,0},{-1,-1},{0,0},{0,-1}});
    mluseAdv_->AssignWENOStencils("edge","(3,3)",{{0,-1},{-2,-1},{-1,0},{-1,-2}});
    mluseAdv_->AssignLinearWgts("edge","(2,2)",{1,1,1,1});
    mluseAdv_->AssignLinearWgts("edge","(3,3)",{5,5,5,5});
    mluseAdv_->UpdateNonLinearWgts(mi, "edge", edge, "HD"); 
    mluseAdv_->UpdateNonLinearWgts(mi, "edge", edge, "CD"); 

    mluseAdv_->AddMLWENOLevel("corner", {"(3,3)", "(2,2)"}, mlpPtr_);
    mluseAdv_->AssignWENOStencils("corner","(2,2)",{{-1,0},{-1,-1},{0,0},{0,-1}});
    mluseAdv_->AssignWENOStencils("corner","(3,3)",{{-2,-2},{0,0},{-2,0},{0,-2}});
    mluseAdv_->AssignLinearWgts("corner","(2,2)",{1,1,1,1});
    mluseAdv_->AssignLinearWgts("corner","(3,3)",{5,5,5,5});
    mluseAdv_->UpdateNonLinearWgts(mi, "corner", corner, "HD"); 
    mluseAdv_->UpdateNonLinearWgts(mi, "corner", corner, "CD"); 

    // --- save diffusion mluse for later ^_^

    //cout <<mluseAdv_->Evaluate({0,-0.5}, {2,2}, mi, "interior", "HD") <<endl;

    return 1;
}

int Driver::UpdateSmoothnessIndicator(){

    mlpPtr_->UpdateSmoothnessIndic(mi, mi.localCD, "CD");
    mlpPtr_->UpdateSmoothnessIndic(mi, mi.localCD, "HD");

    // Diffusion later

    return 1;
}

int Driver::UpdateNonlinearWgts(){

    mluseAdv_->UpdateNonLinearWgts(mi, posSet_, locFuncSet_, fieldSet_);

    return 1;
}

inline vertex effectVel(const vertex& vf, const vertex& vs, const double& c, const double& phif){

    return vf*phif*c + (1-c)*vs;
}

int Driver::SingleEdgeFlux(const indice& localedge,
                           extractEdgeInfoFunc edgeinfo,
                           double& workHD,
                           double& workCD,
                           fluxFunc      fluxfuncAdv, 
                           fluxFuncBndry fluxfuncbndryAdv,
                           fluxFunc      fluxfuncDif,
                           fluxFuncBndry fluxfuncbndryDif){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Compute flux on one given edge
    indice      gCellIn, gCellOut, gCellInside;
    std::string locationIn, locationOut, locationInside;
    edgeEnds<vertex> edgeEndsVertex;
    edgeEnds<indice> edgeEndsIndice;

    // returns values indicating horizontal or vertical edges
    // 1 for horizontal edge, and 2 for vertical edge
    int edgeFlag = edgeinfo(mi, localedge, mi.ghostShiftVertex, 
                            gCellOut, gCellIn, edgeEndsVertex, edgeEndsIndice);

    gCellInside = PickCellInside(mi, gCellIn, gCellOut);

    std::vector<vertex> edge {edgeEndsVertex.start, edgeEndsVertex.end};

    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    // Extract gauss points
    std::vector<vertex> gauss_p;
    gauss_p.resize(gpe.size());

    for (int g=0; g<gpe.size(); g++){
    gauss_p.at(g) = GaussMapPointsEdge({gpe[g]}, edge);}

    // Extract velocity on the given gauss points 
    vector<vertex> vel_darcy = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gauss_p, gCellInside,*hdiv_,*basis_);

    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gauss_p, gCellInside,*br_,*basis_);

    // Computable in and out cell index
    indice comp_gCellIn, comp_gCellOut;
    comp_gCellIn = gCellIn; comp_gCellOut = gCellOut;

    vector<double> HDout, CDout, HDin, CDin;

    HDout.resize(gpe.size());
    CDout.resize(gpe.size());
    HDin.resize(gpe.size());
    CDin.resize(gpe.size());

    // Extract phase behavior on the given gauss points
    vector<double> phifOut, phifIn, cfOut, cfIn;

    phifOut.resize(gpe.size());
    phifIn.resize(gpe.size());
    cfIn.resize(gpe.size());
    cfOut.resize(gpe.size());

    // Compute effective velocity with kappa values
    vector<vertex> velEffectHD, velEffectCD;

    velEffectHD.resize(gpe.size());
    velEffectCD.resize(gpe.size());

    // Equivelent of setting free flow boundary condition
    if (OutBndryCell(mi, gCellIn)){
        // gCellIn out of the boundary
        comp_gCellIn = gCellOut;
    } else if (OutBndryCell(mi, gCellOut)){
        // gCellOut out of the boundary
        comp_gCellOut = gCellIn;
    }

    for (int g=0; g<gpe.size(); g++){
        HDout.at(g) = mluseAdv_->Evaluate(
        gauss_p.at(g), comp_gCellOut, mi, location(mi,comp_gCellOut),"HD");
        CDout.at(g) = mluseAdv_->Evaluate(
        gauss_p.at(g), comp_gCellOut, mi, location(mi,comp_gCellOut),"CD");

        HDin.at(g) = mluseAdv_->Evaluate(
        gauss_p.at(g), comp_gCellIn, mi, location(mi,comp_gCellIn),"HD");
        CDin.at(g) = mluseAdv_->Evaluate(
        gauss_p.at(g), comp_gCellIn, mi, location(mi,comp_gCellIn),"CD");

        // temperatury pressure value
        double depth  = myPhase->pPtr->GetDepth(gauss_p.at(g)[1], myPhase->pp->l0);  
        double lithoP = myPhase->pPtr->GetScaledLithoP(depth);

        // Evaluate phase behavior at outside of the edge
        myPhase->pPtr->evalPhase(HDout.at(g), CDout.at(g), lithoP);

        phifOut.at(g) = myPhase->pPtr->phi.mlt;

        // Compute kappa outside
        double kappaOut = phifOut.at(g) * myPhase->pPtr->Getcf(CDout.at(g)) / CDout.at(g);

        double lambdaOut = phifOut.at(g) * myPhase->pPtr->Getef() / HDout.at(g); 

        // Evaluate phase behaviro at inside of the edge
        myPhase->pPtr->evalPhase(HDin.at(g), CDin.at(g), lithoP);

        phifIn.at(g) = myPhase->pPtr->phi.mlt;

        // Compute kappa outside
        double kappaIn = phifIn.at(g) * myPhase->pPtr->Getcf(CDin.at(g)) / CDin.at(g);

        double lambdaIn = phifOut.at(g) * myPhase->pPtr->Getef() / HDin.at(g);

        double kappa_mean  = harmonic_mean(kappaIn, kappaOut);
        double lambda_mean = harmonic_mean(lambdaIn, lambdaOut);

        double phif_mean = harmonic_mean(phifIn.at(g), phifOut.at(g));

        // Compute effective velocity
        velEffectHD.at(g) = effectVel(vel_darcy.at(g), vel_stokes.at(g), lambda_mean, phif_mean);  
        velEffectCD.at(g) = effectVel(vel_darcy.at(g), vel_stokes.at(g), kappa_mean, phif_mean);  
    }

    // Classify between different flux situation
    if (OutBndryCell(mi, gCellIn)){

        workHD = fluxfuncbndryAdv(mi,gwe, velEffectHD, HDout, unitNormal, len, comp_gCellIn, 
        edgeFlag, "HD");

        workCD = fluxfuncbndryAdv(mi,gwe, velEffectCD, CDout, unitNormal, len, comp_gCellIn, 
        edgeFlag, "CD");

    } else if (OutBndryCell(mi, gCellOut)){

        workHD = fluxfuncbndryAdv(mi,gwe, velEffectHD, HDin, unitNormal, len, comp_gCellOut,
        edgeFlag, "HD");

        workCD = fluxfuncbndryAdv(mi,gwe, velEffectCD, CDin, unitNormal, len, comp_gCellOut,
        edgeFlag, "CD");

    } else {

        workHD = fluxfuncAdv(gwe, velEffectHD, HDin, HDout, unitNormal, len);

        workCD = fluxfuncAdv(gwe, velEffectCD, CDin, CDout, unitNormal, len);

    }

    return 1;
}

int Driver::UpdateEdgeFluxAll(vector<double>& edgefluxHD,
                              vector<double>& edgefluxCD,
                              fluxFunc      fluxfuncAdv, 
                              fluxFuncBndry fluxfuncbndryAdv,
                              fluxFunc      fluxfuncDif,
                              fluxFuncBndry fluxfuncbndryDif){

    // resize two vectors
    edgefluxHD.resize(mi.MPIlocalVertEdgeSize + mi.MPIlocalHoriEdgeSize);
    edgefluxCD.resize(mi.MPIlocalVertEdgeSize + mi.MPIlocalHoriEdgeSize);

    // Assert edge flux vector sizes are correct
    assert(edgefluxHD.size() == mi.MPIlocalVertEdgeSize + mi.MPIlocalHoriEdgeSize);
    assert(edgefluxCD.size() == mi.MPIlocalVertEdgeSize + mi.MPIlocalHoriEdgeSize);

    double fluxHD, fluxCD;

    for (int j=0; j<mi.MPIlocalVertexSize[1]; j++){
        for (int i=0; i<mi.MPIlocalCellSize[0]; i++){
            // Vertical edge computed first
            SingleEdgeFlux({i,j}, extractVertEdgeInfo, fluxCD, fluxHD, fluxfuncAdv, fluxfuncbndryAdv, fluxfuncAdv, fluxfuncbndryAdv);
        }
    }

    for (int j=0; j< mi.MPIlocalCellSize[1]; j++){
        for (int i=0; i<mi.MPIlocalVertexSize[0]; i++){
            // Horizontal edge computed second
            SingleEdgeFlux({i,j}, extractHoriEdgeInfo, fluxCD, fluxHD, fluxfuncAdv, fluxfuncbndryAdv, fluxfuncAdv, fluxfuncbndryAdv);
        }
    }

    return 1;
}

int Driver::ComputeCellFlux(const indice& lCell,
                            double& fluxHD,
                            double& fluxCD,
                            const vector<double>& edgeFluxHD,
                            const vector<double>& edgeFluxCD){

    // Can be modified later if multiple fields are added
    indice left   = lCell;
    indice right  = {lCell[0]+1, lCell[1]};
    indice bottom = lCell;
    indice top    = {lCell[0], lCell[1]+1};

    int left_flat  = FlatIndic(mi.MPIlocalVertexSize[0],left);
    int right_flat = FlatIndic(mi.MPIlocalVertexSize[0],left);
    int bottom_flat= FlatIndic(mi.MPIlocalCellSize[0], bottom) + mi.MPIlocalHoriEdgeSize;
    int top_flat   = FlatIndic(mi.MPIlocalCellSize[0], top) + mi.MPIlocalHoriEdgeSize;       

    fluxHD = edgeFluxHD.at(left_flat) - edgeFluxHD.at(right_flat) + edgeFluxHD.at(bottom_flat) - edgeFluxHD.at(top_flat); 
    fluxHD = edgeFluxCD.at(left_flat) - edgeFluxCD.at(right_flat) + edgeFluxCD.at(bottom_flat) - edgeFluxCD.at(top_flat); 

    return 1;
}
