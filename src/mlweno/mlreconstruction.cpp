// ==========================================================================================
// Class of multi level reconstructions, managing information related to smoothness indicator
// and nonlinear weights between reconstruction levels.
// ==========================================================================================
#include "mlreconstruction.h"

using namespace MLWENO;

//! Add a single weno reconstruction level to the computation
void multiLevelReconstruction::AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY){
    singleLevelReconstruction * slrPtr = new singleLevelReconstruction(stencilSizeX,stencilSizeY);
    slrPtr->CreateStencilPolynomials(mi);
    allLevels_.push_back(slrPtr); 
}

void multiLevelReconstruction::AddWgts_() {

    map<int, double> lw;
    map<int, int> bias;

    const int sizeX = allLevels_[baseReconstMethod_.size()-1]->GetSizeX();

    for (auto& b: baseReconstMethod_[baseReconstMethod_.size()-1]){
             
        lw.insert({FlatIndic(sizeX,b),1.0});
        bias.insert({FlatIndic(sizeX,b),0});
    }

    linearWgts_.push_back(lw);
    etaBias_.push_back(bias);
}

void multiLevelReconstruction::ResetWgts_() {
    // Initialize linear weights and non linear weights
    // with the given information on reconstruction method
    linearWgts_.clear();
    etaBias_.clear();

    linearWgts_.resize(baseReconstMethod_.size());
    etaBias_.resize(baseReconstMethod_.size());

    for (int i=0; i<baseReconstMethod_.size();i++){
        const int sizeX = allLevels_[i]->GetSizeX();
        for (auto& b : baseReconstMethod_[i]){
            linearWgts_[i].insert({FlatIndic(sizeX, b),1.0});
            etaBias_[i].insert({FlatIndic(sizeX,b),0});
        }
    }

}

void multiLevelReconstruction::UpdateNonLinearWgts_(const MeshInfo& mi, indice start){

    assert(baseReconstMethod_.size() == allLevels_.size());

    vector< map<int,double> > nlw(linearWgts_.size());

    double sum = 0.0;

    for (int l=0; l<allLevels_.size(); l++){
       const int sizeX = allLevels_[l]->GetSizeX();
       const int sizeY = allLevels_[l]->GetSizeY();
       for (auto& i:baseReconstMethod_[l]){
            indice owner = start + i;
            if (allLevels_[l]->CheckExist(mi, owner)){
                double scale = allLevels_[l]->GetScale(FlatIndic(mi,owner));
                double sm = allLevels_[l]->CalculateSmoothnessIndic(mi,owner);
                // Get updated smoothness indicators
                double value = linearWgts_[l].at(FlatIndic(sizeX,i))/ 
                               //pow(sm + scale*scale*eps0_ , max(sizeX, sizeY)) * 
                               pow(sm + scale*scale*eps0_ , sizeX+sizeY) * 
                               pow(eps0_*scale / sm+eps0_*
                               scale, etaBias_[l].at(FlatIndic(sizeX,i)));
                nlw[l].insert({FlatIndic(sizeX,i) , value});
                sum += value;
            }            

        }
    }

    for (int l=0; l<allLevels_.size(); l++){
        if (nlw[l].empty() ==0){
            for (auto & in:nlw[l]){
                in.second = in.second/sum; 
            }
        }
    }

    nonLinearWgts_.erase(FlatIndic(mi,start));
    nonLinearWgts_.insert({FlatIndic(mi,start) , nlw});

}

void multiLevelReconstruction::UpdateFirstStageNonLinearWgts_(const MeshInfo& mi, indice start){

    assert(baseReconstMethod_.size() == allLevels_.size());

    vector< map<int,double> > nlw(linearWgts_.size());

    double sum = 0.0;

    for (int l=0; l<allLevels_.size(); l++){
       const int sizeX = allLevels_[l]->GetSizeX();
       const int sizeY = allLevels_[l]->GetSizeY();
       for (auto& i:baseReconstMethod_[l]){
            indice owner = start + i;
            if (allLevels_[l]->CheckExist(mi, owner)){
                double scale = allLevels_[l]->GetScale(FlatIndic(mi,owner));
                double sm = allLevels_[l]->CalculateSmoothnessIndic(mi,owner);
                // Get updated smoothness indicators
                int power = 0;
                if (allLevels_[l]->GetSizeX() * allLevels_[l]->GetSizeY() == 1){
                    power = 1;} else {
                    power = 2;
                }

                double value = linearWgts_[l].at(FlatIndic(sizeX,i))/ 
                               pow(sm + scale*scale*eps0_ , power); 
                nlw[l].insert({FlatIndic(sizeX,i) , value});
                sum += value;
            }            

        }
    }

    for (int l=0; l<allLevels_.size(); l++){
        if (nlw[l].empty() ==0){
            for (auto & in:nlw[l]){
                in.second = in.second/sum; 
            }
        }
    }

    nonLinearWgts_.erase(FlatIndic(mi,start));
    nonLinearWgts_.insert({FlatIndic(mi,start) , nlw});

}

void multiLevelReconstruction::UpdateTwoStageNonLinearWgts_(const MeshInfo& mi, indice start){

    UpdateFirstStageNonLinearWgts_(mi,start);

    vector< map<int,double> > nlw(linearWgts_.size());

    double sum = 0.0;

    for (int l=0; l<allLevels_.size(); l++){
       const int sizeX = allLevels_[l]->GetSizeX();
       const int sizeY = allLevels_[l]->GetSizeY();
       for (auto& i:baseReconstMethod_[l]){
            indice owner = start + i;
            if (allLevels_[l]->CheckExist(mi, owner)){
                double scale = allLevels_[l]->GetScale(FlatIndic(mi,owner));
                double sm = allLevels_[l]->CalculateSmoothnessIndic(mi,owner);
                // Get updated smoothness indicators
                vector< map<int, double> > oldnlw = nonLinearWgts_.at(FlatIndic(mi,start));
                double value = oldnlw[l].at(FlatIndic(sizeX,i))/ 
                               //pow(sm + scale*scale*eps0_ , max(sizeX, sizeY)) * 
                               pow(sm + scale*scale*eps0_ , sizeX+sizeY) * 
                               pow(eps0_*scale / sm+eps0_*
                               scale, etaBias_[l].at(FlatIndic(sizeX,i)));
                nlw[l].insert({FlatIndic(sizeX,i) , value});
                sum += value;
            }            

        }
    }

    for (int l=0; l<allLevels_.size(); l++){
        if (nlw[l].empty() ==0){
            for (auto & in:nlw[l]){
                in.second = in.second/sum; 
            }
        }
    }

    nonLinearWgts_.erase(FlatIndic(mi,start));
    nonLinearWgts_.insert({FlatIndic(mi,start) , nlw});
}

void multiLevelReconstruction::UpdateNonLinearWgts(const MeshInfo& mi){
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){
        indice add {i,j};
        indice start = mi.MPIlocalCellStart+add;
        UpdateNonLinearWgts_(mi,start);
    } }
}

void multiLevelReconstruction::UpdateTwoStageNonLinearWgts(const MeshInfo& mi){
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){
        indice add {i,j};
        indice start = mi.MPIlocalCellStart+add;
        UpdateTwoStageNonLinearWgts_(mi,start);
    } }
}

void multiLevelReconstruction::GetInfo(){

    cout << "There are " <<allLevels_.size()<< " levels." << endl;

    for (int l=0; l<allLevels_.size(); l++){
        allLevels_[l]->CheckStencils(); 
    }
}

// ===================================================================================
void multiLevelReconstruction::PrintSmoothnessIndicator(const MeshInfo& mi){

    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){

        indice add {i,j};
        indice start = mi.MPIlocalCellStart + add;

        cout << "Reconstruction at cell ( " << start[0] << ", " << start[1] << ")" << endl; 
        for (int l=0; l < allLevels_.size(); l++) {
            allLevels_[l]->PrintSmoothnessIndicator(mi);
        }
 
    }}

}

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
