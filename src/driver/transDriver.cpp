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

int Driver::SingleEdgeFlux_(const indice& localedge,
                            extractEdgeInfoFunc edgeinfo,
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
    CDin.resize(gep.size());

    // Extract phase behavior on the given gauss points
    vector<double> phifOut, phifIn, cfOut, cfIn, kappaIn, kappaOut;

    phifOut.resize(gpe.size());
    phifIn.resize(gpe.size());
    cfIn.resize(gpe.size());
    cfOut.resize(gpe.size());

    // Compute effective velocity with kappa values
    vector<vertex> velEffectOut, velErrectIn;

     
    // Equivelent of setting free flow boundary condition
    if (OutBndryCell(mi, globalCellIn)){
        // gCellIn out of the boundary
        comp_gCellIn = gCellOut;
    } else if (OutBndryCell(mi, globalCellOut)){
        // gCellOut out of the boundary
        comp_gCellOut = gCellIn;
    }

    for (int g=0; g<gpe.size(); g++){
        HDout.at(g) = mluseAdv_.Evaluate(
        gauss_p.at(g), comp_gCellOut, mi, location(mi,comp_gCellOut),"HD");
        CDout.at(g) = mluseAdv_.Evaluate(
        gauss_p.at(g), comp_gCellOut, mi, location(mi,comp_gCellOut),"CD");

        HDin.at(g) = mluseAdv_.Evaluate(
        gauss_p.at(g), comp_gCellIn, mi, location(mi,comp_gCellIn),"HD");
        CDin.at(g) = mluseAdv_.Evaluate(
        gauss_p.at(g), comp_gCellIn, mi, location(mi,comp_gCellIn),"CD");

        double depth  = phase->pPtr->GetDepth(gauss_p.at(g)[1]);  
        double lithoP = phase->pPtr->GetScaledLithoP(depth);

        // Evaluate phase behavior at outside of the edge
        phase->pPtr->evalPhase(HDout.at(g), CDout.at(g), lithoP);

        phifOut.at(g) = phase->pPtr->phi.mlt;
        cfOut.at(g)   = CDOut.at(g) - phase->pPtr->phi.opx;

        // Compute kappa outside
        kappaOut.at(g) = phiOut.at(g) * cfOut.at(g) / CDOut.at(g);

        // Evaluate phase behaviro at inside of the edge
        phase->pPtr->evalPhase(HDin.at(g), CDin.at(g), lithoP);

        phifIn.at(g) = phase->pPtr->phi.mlt;
        cfIn.at(g)   = CDin.at(g) - phase->pPtr->phi.opx;

        // Compute kappa outside
        kappaIn.at(g) = phiIn.at(g) * cfIn.at(g) / CDIn.at(g);
    }

   

    return 1;
}
