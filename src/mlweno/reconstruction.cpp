#include "reconstruction.h"

void reconstruction::UpdateStencilSizeMax_(int[2] newStencilSize){

    if (newStencilSize[0]*newStencilSize[1] > stencilSizeMax[0]*stencilSizeMax[1]){
        stencilSizeMax_[0] = newStencilSize[0];
        stencilSizeMax_[1] = newStencilSize[1];
    }

}

void reconstruction::AddStencil(int[2] stencilSize, int[2] shift) {

    stencilSize_.push_back(stencilSize);
    shift_.push_back(shift);

    stencilNum_ = stencilSize_.size();

    UpdateStencilSizeMax_(stencilSize);

    // Create stencil Indice according to the information
    stencil <indice> siNow;

     

    stencilIndice_.push_back(siNow);
}


