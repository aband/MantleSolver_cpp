#include "reconstruction.h"

using namespace MLWENO;

void reconstruction::UpdateStencilSizeMax_(int* newStencilSize){

    if (newStencilSize[0]*newStencilSize[1] > stencilSizeMax_[0]*stencilSizeMax_[1]){
        stencilSizeMax_[0] = newStencilSize[0];
        stencilSizeMax_[1] = newStencilSize[1];
    }

}

void reconstruction::AddStencil(int* stencilSize, indice shift) {

    stencilSize_.push_back(stencilSize);
    shift_.push_back(shift);

    stencilNum_ = stencilSize_.size();

    UpdateStencilSizeMax_(stencilSize);

    // Create stencil Indice according to the information
    stencil <indice> siNow(stencilSize[0],stencilSize[1]);

    for (int j=0; j<stencilSize[1]; j++){
    for (int i=0; i<stencilSize[0]; i++){

        indice cell = {i,j};

        siNow(i,j) = cell + shift;

    }}

    stencilIndice_.push_back(siNow);
}

void reconstruction::AddStencil(vector<int*> stencilSizes, vector<indice> shifts) {

    assert(stencilSizes.size() == shifts.size());

    for (int i=0; i<stencilSizes.size(); i++){
        AddStencil(stencilSizes[i], shifts[i]); 
    }

}

void reconstruction::PrintStencils() const {

    for (int s=0; s<stencilIndice_.size(); s++){
        stencil <indice> siNow = stencilIndice_[s];
        for (int j=0; j<siNow.getJ(); j++){
        for (int i=0; i<siNow.getI(); i++){
            indice iNow = siNow(i,j);
            cout << std::setw(2) << "(" << iNow[0] << "," << iNow[1] << ") " ;
        } cout << endl;}
        cout << endl;
    }

}

void reconstruction::Clear() {

    shift_.clear(); 

    stencilSizeMax_[0] = 0;
    stencilSizeMax_[1] = 0;

    stencilNum_ = 0;

    stencilIndice_.clear();
    stencilPolyn_.clear();

}
