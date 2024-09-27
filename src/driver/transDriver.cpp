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
