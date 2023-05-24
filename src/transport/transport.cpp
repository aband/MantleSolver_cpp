#include "transport.h"

// ========== Transport ===========================================
// Inherit from advection, diffusion and reaction class

void transport::AddLevel(const MeshInfo& mi, const int& stencilSizeX,
                                             const int& stencilSizeY){
    mlpPtr_->AddLevel(mi, stencilSizeX, stencilSizeY);
}

void transport::UpdateSmoothnessIndic(const MeshInfo& mi){
    mlpPtr_->UpdateSmoothnessIndic(mi);    
}

void transport::UpdateSmoothnessIndicAndDerivative(const MeshInfo& mi){
    mlpPtr_->UpdateSmoothnessIndic(mi); 
    mlpPtr_->UpdateDerivSmoothnessIndic(mi);
}

void transport::AssignReconstMethod(unordered_map<std::string, vector<indice>>& reconstMethods,
                                    std::string key, const vector<indice>& brm){
    //! Make sure the reconstruction key is new
    assert(reconstMethods.count(key) == 0);
    //! Add this new pair to reconstMethods
    reconstMethods.insert(std::pair<std::string, vector<indice>>(key, brm));
}

void transport::CreateMLWENO(const MeshInfo& mi){
    advection::CreateMLWENO((*mlpPtr_), mi);
    diffusion::CreateMLWENO((*mlpPtr_), mi);
}
