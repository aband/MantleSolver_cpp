// ==========================================================================================
// Class of multi level reconstructions, managing information related to smoothness indicator
// and nonlinear weights between reconstruction levels.
// ==========================================================================================
#include "mlreconstruction.h"

using namespace MLWENO;

//! Add a single weno reconstruction level to the computation
void multiLevelReconstruction::AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY,
                                        vector<indice> brm){
    //! Create single level key 
    std::string key = '(' + std::to_string(stencilSizeX) + ',' + 
                            std::to_string(stencilSizeY) + ')';

    //! Check if the new level has never been defined before
    assert(reconstMethods_.count(key) == 0);

    //! Insert new key into allLevels;
    wenoLevels_.insert(key);

    //! Initialize single level reconstruction class pointer.
    singleLevelReconstruction * slrPtr = new singleLevelReconstruction(stencilSizeX,stencilSizeY);
    slrPtr->CreateStencilPolynomials(mi);

    //! Create map from single level key to created single level reconstrucion
    reconstLevels_.insert(std::pair<std::string,singleLevelReconstruction *>(key,slrPtr));

    //! Create map from single level key to reconstruction method (how to find stencils)
    reconstMethods_.insert(std::pair<std::string, vector<indice>>(key,brm)); 

    //! Update lowest and highest order of reconstruction level
    lowestLevel_ = reconstMethods_.begin()->first;
    highestLevel_ = reconstMethods_.rbegin()->first;

    //! Accumulate total levels
    totalLevels_ += brm.size();
}

//! Specify boundary layers (cells near boundary that need additional reconstruciton level then interior cells)
void multiLevelReconstruction::SeparateBoundaryLayer(const MeshInfo& mi, const int& layerSize){
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){
        indice shift {i,j}; 
        indice global = mi.MPIlocalCellStart + shift;
        if (global[0] == 0 + layerSize-1 || global[0] == mi.MPIglobalCellSize[0] - layerSize ||
            global[1] == 0 + layerSize-1 || global[1] == mi.MPIglobalCellSize[1] - layerSize){
            boundaryCells_.insert(FlatIndic(mi,global));
        } else {
            interiorCells_.insert(FlatIndic(mi,global));
        }
    }}

    //! Create two different weno reconstruction levels for interior and boundary cells.
    SeparateReconstMethods_();
}

void multiLevelReconstruction::SeparateBoundaryLayer(const MeshInfo& mi){
    SeparateBoundaryLayer(mi,1);
}

void multiLevelReconstruction::SeparateReconstMethods_(){
    if (lowestLevel_ != "(1,1)"){
        cout << "Lowest reconstruction level is not constant." << endl;
        interiorLevels_ = wenoLevels_;
        boundaryLevels_ = wenoLevels_;
    } else {
        interiorLevels_ = wenoLevels_;
        boundaryLevels_ = wenoLevels_;
        interiorLevels_.erase(lowestLevel_);
    }

}

//! Update non linear weights for all cell reconstructions.
void multiLevelReconstruction::UpdateOneStageNonLinearWgts_(const MeshInfo& mi, int flatGlobal,
                                                            unordered_set<std::string> levels){

    unordered_map<std::string, unordered_map<int, double>> nlw;

    double sum = 0.0;

    for (auto const& level : levels){
        const int sizeX = reconstMethods_[level]->getsizeX();
        const int sizeY = reconstMethods_[level]->getsizeY();
        for (auto const& rm : reconstMethods_[level]){
            indice owner = Bend(mi,flatGlobal) + rm;
            if (reconstLevels_[rm]->CheckExist(mi,owner)){
                double scale = reconstLevels_[level]->GetScale(FlatIndic(mi, owner));
                double sm = reconstLevels_[level]->GetSmoothnessIndic(mi, owner);
                double value = 1.0/pow(sm + scale*scale*eps0_ , sizeX+sizeY) * 
                               pow(eps0_*scale / sm+eps0_*
                               scale, etaBias_[l].at(FlatIndic(sizeX,i)));
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                nlw[] 
                sum += value;
            }
        }
    }

}

void multiLevelReconstruction::UpdateOneStageNonLinearWgts(const MeshInfo& mi){

    for (auto const& singleLevel : wenoLevels_){
        reconstLevels_[singleLevel]->UpdateSmoothnessIndic(mi);
    }

    for (auto const& it: interiorCells_){
        UpdateOneStageNonLinearWgts_(mi, it, interiorLevels_);
    }

    for (auto const& it: boundaryCells_){
        UpdateOneStageNonLinearWgts_(mi, it, boundaryLevels_);
    }

}

void multiLevelReconstruction::UpdateTwoStageNonLinearWgts(const MeshInfo& mi){


}

// ===================================================================================
void multiLevelReconstruction::PrintBoundaryLayer(const MeshInfo& mi){
    for (auto const& it : interiorCells_){
        indice global = Bend(mi,it); 
        cout << "(" << global[0] << "," << global[1] << ") "; 
    }cout << endl;

    for (auto const& it : boundaryCells_){
        indice global = Bend(mi,it); 
        cout << "(" << global[0] << "," << global[1] << ") "; 
    }cout << endl;

}


void multiLevelReconstruction::GetInfo(){
    //! Print added levels and reconstruction methods
    cout << "There are " <<reconstLevels_.size()<< " levels pre computed." << endl;

    for (auto const& it : reconstLevels_){
        cout << it.first << " " ; 
        (it.second)->CheckStencils();
    }

    //cout << "Highest order of level is " << highestLevel_ << endl;
    //cout << "Lowest order of level is " << lowestLevel_ << endl;

}

void multiLevelReconstruction::PrintSmoothnessIndicator(const MeshInfo& mi){


        for (auto const& singleLevel : wenoLevels_) {
            cout << "Current reconstruction level is :" << singleLevel << endl;
            reconstLevels_[singleLevel]->PrintSmoothnessIndicator(mi);
        }
 
/*
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){

        indice add {i,j};
        indice start = mi.MPIlocalCellStart + add;

        cout << "Reconstruction at cell ( " << start[0] << ", " << start[1] << ")" << endl; 
        for (auto const& singleLevel : wenoLevels_) {
            cout << "Current reconstruction level is :" << singleLevel << endl;
            reconstLevels_[singleLevel]->PrintSmoothnessIndicator(mi);
        }
 
    }}
*/

}


/*

void multiLevelReconstruction::PrintNonLinearWgts(const MeshInfo& mi){

    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){

        indice add {i,j};
        indice start = mi.MPIlocalCellStart + add;
        const vector<map<int,double>>& nlw = nonLinearWgts_[FlatIndic(mi,start)];

        cout << "Reconstruction at cell ( " << start[0] << ", " << start[1] << ")" << endl; 
        for (int l=0; l < allLevels_.size(); l++) {
            int sizeX = allLevels_[l]->GetSizeX();

            if (nlw[l].empty() == 0) {
                for (auto & in:nlw[l]){
                    indice m = Bend(sizeX,in.first);
                    cout << "At Level " << l << " reconstruction at ( " << m[0] << ", "
                         << m[1] << ") " << " with wgt " << in.second << endl;
                }
            }
        }
    }}
}

void multiLevelReconstruction::Clear(){
    for (int i=0; i<allLevels_.size(); i++){
        delete allLevels_[i];
    }
    allLevels_.clear();
}
*/
