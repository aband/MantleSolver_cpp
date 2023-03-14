#include "transport.h"

advection::SelectReconstLevels(unordered_set<std::string> reconstLevels){

    
}

advection::Func_(double x, double y, double u, double t){
    return func(double x, double y, double u, double t);
}

advection::dFunc_(double x, double y, double u, double t){
    return dfunc(double x, double y, double u, double t);
}


transport::AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY, vector<indice> brm){
    mlrPtr_->AddLevel(mi,stencilSizeX,stencilSizeY,brm);
}


