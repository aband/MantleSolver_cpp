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

int Driver::PrepareDefaultMLWENO(){

    // Prepare levels 
    AddLevels(1);
    AddLevels(2);
	 AddLevels(3);

    AddLevels(4,3);
    AddLevels(3,4);

    // Create Smoothness Indicator for the first time
    // mlpPtr_->UpdateSmoothnessIndic(mi, mi.localCD, "CD");
    mlpPtr_->UpdateSmoothnessIndic(mi, mi.localHD, "HD");

    mluseAdv_->AddMLWENOLevel("interior", {"(3,3)","(2,2)"}, mlpPtr_);

    mluseAdv_->AssignWENOStencils("interior","(2,2)",{{-1,0},{-1,-1},{0,0},{0,-1}});
    mluseAdv_->AssignWENOStencils(0,"(3,3)",{{-1,-1}});

    mluseAdv_->AssignLinearWgts("interior","(2,2)",{1,1,1,1});
    mluseAdv_->AssignLinearWgts("interior","(3,3)",{10});

    mluseAdv_->UpdateNonLinearWgts(mi, "interior", interior, "HD"); 
    mluseAdv_->UpdateNonLinearWgts(mi, "interior", interior, "CD"); 

    // --- save diffusion mluse for later ^_^

    return 1;
}
