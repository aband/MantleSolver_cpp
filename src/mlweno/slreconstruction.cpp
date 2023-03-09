#include "slreconstruction.h"
//! Routines for singleLevelReconstruction

using namespace MLWENO;

singleLevelReconstruction::singleLevelReconstruction(int stencilSizeX, int stencilSizeY){
    stencilSizeX_ = stencilSizeX;
    stencilSizeY_ = stencilSizeY;
 
    //! Create stencil of indice used in creating stencil polynomials
    stencilIndice_.SetStencil(stencilSizeX_, stencilSizeY_); 

    for (int j=0; j<stencilSizeY_; j++){
    for (int i=0; i<stencilSizeX_; i++){
        stencilIndice_(i,j) = {i,j};
    }}
}

//! Calculate smoothness indicator for the single stencil
double singleLevelReconstruction::CalculateSmoothnessIndic(const MeshInfo& mi, indice owner) {
    return singleLevel_[FlatIndic(mi, owner)]->GetSmoothIndic(mi,stencilIndice_);
}

void singleLevelReconstruction::IdentifyInteriorCell_(const MeshInfo& mi){
    //! Should be called each time when a new level is added to reconstruction
    for (int j=-2; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=-2; i<mi.MPIlocalCellSize[0]; i++){ 
       if (j+mi.MPIlocalCellStart[1]+stencilSizeY_-1<mi.MPIglobalCellSize[1] &&
           i+mi.MPIlocalCellStart[0]+stencilSizeX_-1<mi.MPIglobalCellSize[0] &&
           j+mi.MPIlocalCellStart[1]>-1 &&
           i+mi.MPIlocalCellStart[0]>-1 ){

               interior_.insert(FlatIndic(mi,i,j));
           }
    }}
}

/**!
 * Two steps of creating stencil polynomials for a given reconstruction level.
 *     Step 1: Identify all the stencils covering entire computational domain.
 *     Step 2: Compute stencil polynomials for each stencil.
 */
void singleLevelReconstruction::CreateStencilPolynomials(const MeshInfo& mi){
    IdentifyInteriorCell_(mi);
    ComputeStencilPolyn_(mi);
}

//! Calculate center vertex of a given stencil
const vertex singleLevelReconstruction::ComputeStencilCenter_(const MeshInfo& mi, int flat){
    vertex work  = {0.0,0.0};
    indice original = Bend(mi,flat);

    original[0] = original[0] + mi.vertexGhostLayerSize;
    original[1] = original[1] + mi.vertexGhostLayerSize;

    work += mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0],original)];

    original[0] = original[0] + stencilSizeX_;
 
    work += mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0],original)];

    original[1] = original[1] + stencilSizeY_;
 
    work += mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0],original)];

    original[0] = original[0] - stencilSizeX_;
 
    work += mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0],original)];

    return work/4.0;
}

//! Update smoothness indicators for the entire level
void singleLevelReconstruction::UpdateSmoothnessIndic(const MeshInfo& mi){
    for (auto& ind:interior_){
        smoothnessIndic_[ind] = singleLevel_[ind]->GetSmoothIndic(mi,stencilIndice_);
    }
}

//! Compute stencil polynomials for the entire level
void singleLevelReconstruction::ComputeStencilPolyn_(const MeshInfo& mi){
    for (auto& flat: interior_){
        singleLevel_[flat] = new stencilPolynomial(Bend(mi,flat), ComputeStencilCenter_(mi,flat));
        singleLevel_[flat]->SetUpScale(mi,stencilIndice_);
        singleLevel_[flat]->SetStencilPolynomials(mi,stencilIndice_);
    }
}

//! Extract smoothness indicator from pre-calculated values
double singleLevelReconstruction::GetSmoothnessIndic(const MeshInfo& mi, indice owner){
    return smoothnessIndic_[FlatIndic(mi,owner)];
}

// ======================================================================
//! Functions checking created stencil polynomials for a single stencil
void singleLevelReconstruction::CheckStencilPolynomials(const MeshInfo& mi, indice start){
    singleLevel_[FlatIndic(mi,start)]->printCoef();
}

//! Print all smoothnessIndicator for all stencils belonging to this level
void singleLevelReconstruction::PrintSmoothnessIndicator(const MeshInfo& mi){

    for (auto& ind:interior_){
        cout << smoothnessIndic_[ind] << endl;
    }

}
