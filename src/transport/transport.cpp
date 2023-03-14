#include "transport.h"

void advection::SelectReconstLevels(unordered_set<std::string> reconstLevels){
    reconstLevels_ = recosntLevels; 
}

//! Return function value defined in the separate function file
void advection::FuncX_(double x, double y, double u, double t){
    return funcX(double x, double y, double u, double t);
}

void advection::dFuncX_(double x, double y, double u, double t){
    return dfuncX(double x, double y, double u, double t);
}

void advection::FuncY_(double x, double y, double u, double t){
    return funcY(double x, double y, double u, double t);
}

void advection::dFuncY_(double x, double y, double u, double t){
    return dfuncY(double x, double y, double u, double t);
}

void diffusion::SelectReconstLevels(unordered_set<std::string> reconstLevels){
    reconstLevels_ = reconstLevels;
}

void reaction::SelectReconstLevels(unordered_set<std::string> reconstLevels){
    reconstLevels_ = reconstLevels;
}

// ==============================================================================================
void transport::AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY, vector<indice> brm){
    mlrPtr_->AddLevel(mi,stencilSizeX,stencilSizeY,brm);
}

void transport::AddReconstMethod(unordered_map<std::string, vector<indice>>& reconstMethods,
                            std::string key, vector<indice> brm){
    //! Make sure the reconstruction key is new
    assert(reconstMethods[key].empty() == 0);
    //! Add this new pair to reconstMethods
    reconstMethods.insert(std::pair<std::string, vector<indice>>(key, brm));

}

unordered_set transport::CreateWenoLevel(const unordered_map<std::string, vector<indice>>& reconstMethods){
    unordered_set work;
    for (auto rm : reconstMethods){
        work.insert(rm.first); 
    }
    return work;
}
