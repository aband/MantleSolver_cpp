#include "reconstMLWENO.h"

using namespace MLWENO;
using namespace tensorProductPoly;

//! A constructor
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

//! A destructor
singleLevelReconstruction::~singleLevelReconstruction(){
    interior_.clear();
    smoothnessIndic_.clear();

    for (auto& it: singleLevel_){
        delete it.second;
    }
    singleLevel_.clear();
}

//! Calculate smoothness indicator for the single stencil
double singleLevelReconstruction::CalculateSmoothnessIndic(const MeshInfo& mi, indice owner) {
    return singleLevel_[FlatIndic(mi, owner)]->GetSmoothIndic(mi,stencilIndice_);
}

void singleLevelReconstruction::IdentifyInteriorCell_(const MeshInfo& mi){

    // Function used to find interior stencils.
    // Has nothing to do with physical domain.

    for (int j=mi.MPIlocalCellStart[1]-mi.cellGhostLayerSize; 
         j<mi.MPIlocalCellStart[1]+mi.MPIlocalCellSize[1]+mi.cellGhostLayerSize; j++){
    for (int i=mi.MPIlocalCellStart[0]-mi.cellGhostLayerSize; 
         i<mi.MPIlocalCellStart[0]+mi.MPIlocalCellSize[0]+mi.cellGhostLayerSize; i++){

        if (j+stencilSizeY_-1 < mi.MPIglobalCellSize[1] // The top of stencil does not exceed maximum Y 
        &&  i+stencilSizeX_-1 < mi.MPIglobalCellSize[0] // The right side of stencil does not exceed maximum X
        &&  j > -1                                      // The bottom of stencil does not exceed minimum Y
        &&  i > -1                                      // The bottom of stencil does not exceed minimum X
        &&  i+stencilSizeX_-1 < mi.MPIlocalCellStart[0]+mi.MPIlocalCellSize[0]+mi.cellGhostLayerSize		  
        // The right side of stencil does not exceed mesh partition
        &&  j+stencilSizeY_-1 < mi.MPIlocalCellStart[1]+mi.MPIlocalCellSize[1]+mi.cellGhostLayerSize		  
        // The top side of stencil does not exceed mesh partition
        ){
            //if(rank == 1){cout << j << " " << i << endl;}
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

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    vertex work  = {0.0,0.0};
    indice original = Bend(mi,flat);

    original[0] = original[0] + mi.vertexGhostLayerSize - mi.MPIlocalCellStart[0];
    original[1] = original[1] + mi.vertexGhostLayerSize - mi.MPIlocalCellStart[1];

//    if (rank == 1){cout << "Compute Stencil Center :" <<mi.lmesh.size() << " " << original[0] << " " << original[1] << " " << 
//				mi.MPIlocalVertexSizeFull[0] << " " << 
//				FlatIndic(mi.MPIlocalVertexSizeFull[0],original) << endl;}

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

void singleLevelReconstruction::UpdateDerivSmoothnessIndic(const MeshInfo& mi){
    for (auto& ind:interior_){
        smoothnessIndicDeriv_[ind] = singleLevel_.at(ind)->GetDerivSmoothIndic(mi,stencilIndice_);
    }
}

//! Compute stencil polynomials for the entire level
void singleLevelReconstruction::ComputeStencilPolyn_(const MeshInfo& mi){
    for (auto& flat: interior_){
        singleLevel_[flat] = new stencilPolynomial(Bend(mi,flat), ComputeStencilCenter_(mi,flat));
        singleLevel_[flat]->SetUpScale(mi,stencilIndice_);
        singleLevel_[flat]->SetStencilPolynomials(mi,stencilIndice_);
        //singleLevel_[flat]->printCoef();
    }
}

/**
 * Evaluate at the given single reconstruction level
 */
double singleLevelReconstruction::Evaluate(const MeshInfo& mi, const indice& owner, const vertex& point){

    if (CheckExist(mi,owner)){
        return singleLevel_[FlatIndic(mi,owner)]->eval(point);
    } else {
        return 0;
    }
}

/**
 * Evaluate at the given single reconstruction level 
 * but evaluate individual stencil polynomials separately not the collapsed one.
 */
double singleLevelReconstruction::Evaluate(const MeshInfo& mi, const indice& owner, const vertex& point, const int& local){

    if (CheckExist(mi,owner)){
        return singleLevel_[FlatIndic(mi,owner)]->eval(point,local);
    } else {
        return 0;
    }
    //return CheckExist(mi, owner) * singleLevel_[FlatIndic(mi,owner)]->eval(point, local);
}

//! Extract smoothness indicator from pre-calculated values
double singleLevelReconstruction::GetSmoothnessIndic(const MeshInfo& mi, indice owner){
    return smoothnessIndic_[FlatIndic(mi,owner)];
}

unordered_map<int,double> singleLevelReconstruction::GetSmoothnessIndicDeriv(const MeshInfo& mi, const indice& owner){
    return smoothnessIndicDeriv_.at(FlatIndic(mi,owner));
}

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

// ============= Multi level reconstruction preparation ========================
MLWENOPrepare::~MLWENOPrepare(){
    for (auto& it : allLevels_){
        delete it.second;
    }
    allLevels_.clear();
}

void MLWENOPrepare::AddLevel(const MeshInfo& mi, const int& stencilSizeX,
                                                 const int& stencilSizeY){

    //! Create single level key 
    std::string key = '(' + std::to_string(stencilSizeX) + ',' + 
                            std::to_string(stencilSizeY) + ')';

    //! Check if the new level has never been defined before
    assert(allLevels_.count(key) == 0);

    //! Initialize single level reconstruction class pointer.
    singleLevelReconstruction * slrPtr = new singleLevelReconstruction(stencilSizeX,stencilSizeY);
    slrPtr->CreateStencilPolynomials(mi);

    //! Create map from single level key to created single level reconstrucion
    allLevels_.insert(std::pair<std::string,singleLevelReconstruction *>(key,slrPtr));

}

void MLWENOPrepare::UpdateSmoothnessIndic(const MeshInfo& mi){

    for (auto const& singleLevel : allLevels_){
        singleLevel.second->UpdateSmoothnessIndic(mi);
        //singleLevel.second->PrintSmoothnessIndicator(mi);
    }
}

void MLWENOPrepare::UpdateDerivSmoothnessIndic(const MeshInfo& mi){
    for (auto const& singleLevel : allLevels_){
        singleLevel.second->UpdateDerivSmoothnessIndic(mi);
    }
}

void MLWENOPrepare::PrintInfo(){

    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    cout << "Current rank is :" << rank << endl;

    //! Print added levels and reconstruction methods
    cout << "There are " <<allLevels_.size()<< " levels created." << endl;

    for (auto const& it : allLevels_){
        cout << it.first << " " ; 
        (it.second)->CheckStencils();
    }
}

// ========= MultiLevelReconstruction ====================================
void multiLevelReconstruction::UpdateNonLinearWgts(const MeshInfo& mi){

}



void multiLevelReconstruction::UpdateNonLinearWgtsSingleCell(const MeshInfo& mi,
                                                             const int& globalCell){


}
