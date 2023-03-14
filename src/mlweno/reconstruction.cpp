#include "reconstruction.h"

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
    //! Should be called each time when a new level is added to reconstruction
    for (int j=-2; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=-2; i<mi.MPIlocalCellSize[0]; i++){
       //! Transform local indices to global indices
       int globali = i+mi.MPIlocalCellStart[0];
       int globalj = j+mi.MPIlocalCellStart[1];

       if (globalj+stencilSizeY_-1<mi.MPIglobalCellSize[1] &&
           globali+stencilSizeX_-1<mi.MPIglobalCellSize[0] &&
           globalj > -1 && globali > -1 ){

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

/**
 * Evaluate at the given single reconstruction level
 */
double singleLevelReconstruction::Evaluate(const MeshInfo& mi, indice owner, vertex point){
    return CheckExist(mi, owner) * singleLevel_[FlatIndic(mi,owner)]->eval(point);
}

//! Extract smoothness indicator from pre-calculated values
double singleLevelReconstruction::GetSmoothnessIndic(const MeshInfo& mi, indice owner){
    return smoothnessIndic_[FlatIndic(mi,owner)];
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

// ============= Multi level reconstruction ====================================
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

    //! Clear old boundary and interior levels
    boundaryLevels_.clear();
    interiorLevels_.clear();

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
                                                            const unordered_set<std::string>& levels){

    unordered_map<std::string, unordered_map<int, double>> nlw;

    double sum = 0.0;

    for (auto const& level : levels){
        const int sizeX = reconstLevels_[level]->GetSizeX();
        const int sizeY = reconstLevels_[level]->GetSizeY();
        for (auto const& rm : reconstMethods_[level]){
            indice owner = Bend(mi,flatGlobal) + rm;
            if (reconstLevels_[level]->CheckExist(mi,owner)){
                double scale = reconstLevels_[level]->GetScale(FlatIndic(mi, owner));
                double sm = reconstLevels_[level]->GetSmoothnessIndic(mi, owner);
                double value = 1.0/pow(sm + scale*scale*eps0_ , sizeX+sizeY);
                nlw[level].insert(std::pair<int, double>(FlatIndic(sizeX,rm), value));
                sum += value;
            }
        }
    }

    for (auto const& level: levels){
        if (nlw[level].empty() == 0){
            for (auto & in: nlw[level]){
                in.second = in.second/sum;
            }
        }
    }

    nonLinearWgts_.erase(flatGlobal);
    nonLinearWgts_.insert(std::pair<int, unordered_map<std::string, unordered_map<int, double>>>(flatGlobal,nlw));

}

void multiLevelReconstruction::UpdateTwoStageNonLinearWgts_(const MeshInfo& mi, int flatGlobal,
                                                            const unordered_set<std::string>& levels){

    unordered_map<std::string, unordered_map<int, double>> tmp;
    unordered_map<std::string, unordered_map<int, double>> nlw;

    double sum = 0.0;

    //! First stage of the calculation of non linear weights
    for (auto const& level : levels){
        const int sizeX = reconstLevels_[level]->GetSizeX();
        const int sizeY = reconstLevels_[level]->GetSizeY();
        for (auto const& rm : reconstMethods_[level]){
            indice owner = Bend(mi,flatGlobal) + rm;
            if (reconstLevels_[level]->CheckExist(mi,owner)){
                double scale = reconstLevels_[level]->GetScale(FlatIndic(mi, owner));
                double sm = reconstLevels_[level]->GetSmoothnessIndic(mi, owner);

                int power = 2;
                double linearWgt = 1.0;
                if (sizeX*sizeY == 1){power = 1; linearWgt = 1e-2;}

                double value = linearWgt/pow(sm + scale*scale*eps0_ , power);
                tmp[level].insert(std::pair<int, double>(FlatIndic(sizeX,rm), value));
                sum += value;
            }
        }
    }

    for (auto const& level: levels){
        if (tmp[level].empty() == 0){
            for (auto & in: tmp[level]){
                in.second = in.second/sum;
            }
        }
    }

    //! Second stage of the calculatino of non linear weights
    sum = 0.0; 
    for (auto const& level : levels){
        const int sizeX = reconstLevels_[level]->GetSizeX();
        const int sizeY = reconstLevels_[level]->GetSizeY();
        unordered_map wgts = tmp[level];

        for (auto const& rm : reconstMethods_[level]){
            indice owner = Bend(mi,flatGlobal) + rm;
            if (reconstLevels_[level]->CheckExist(mi,owner)){
                double scale = reconstLevels_[level]->GetScale(FlatIndic(mi, owner));
                double sm = reconstLevels_[level]->GetSmoothnessIndic(mi, owner);

                double value = wgts[FlatIndic(sizeX,rm)]/pow(sm + scale*scale*eps0_ , sizeX*sizeY);
                nlw[level].insert(std::pair<int, double>(FlatIndic(sizeX,rm), value));
                sum += value;
            }
        }
    }

    for (auto const& level: levels){
        if (nlw[level].empty() == 0){
            for (auto & in: nlw[level]){
                in.second = in.second/sum;
            }
        }
    }

    nonLinearWgts_.erase(flatGlobal);
    nonLinearWgts_.insert(std::pair<int, unordered_map<std::string, unordered_map<int, double>>>(flatGlobal,nlw));

}

void multiLevelReconstruction::UpdateNonLinearWgts(const MeshInfo& mi, const int stage){

    for (auto const& singleLevel : wenoLevels_){
        reconstLevels_[singleLevel]->UpdateSmoothnessIndic(mi);
    }

    if (stage == 1){
        for (auto const& it: interiorCells_){
            UpdateOneStageNonLinearWgts_(mi, it, interiorLevels_);
        }

        for (auto const& it: boundaryCells_){
            UpdateOneStageNonLinearWgts_(mi, it, boundaryLevels_);
        }

    }else if (stage == 2){
        for (auto const& it: interiorCells_){
            UpdateTwoStageNonLinearWgts_(mi, it, interiorLevels_);
        }

        for (auto const& it: boundaryCells_){
            UpdateTwoStageNonLinearWgts_(mi, it, boundaryLevels_);
        }

    }

}

/**
 * Modify reconst method with given key.
 * Only modify existing reconst method.
 * Please use AddLevel to add reconstruction levels.
 */
void multiLevelReconstruction::ModifyReconstMethod(std::string key, vector<indice> newReconstMethod){
    //! Make sure it is a existing key
    assert(reconstMethods_.count(key) != 0); 

    //! Ereasing the existing key
    reconstMethods_.erase(key);

    //! Assign key and reconstruction method pair
    reconstMethods_.insert(std::pair<std::string, vector<indice>>(key, newReconstMethod));
}

/**
 * Modify wenoLevels_.
 * Selecting predefined reconstruction method in the actual computation process.
 */
void multiLevelReconstruction::SelectWenoReconstLevel(unordered_set<std::string> keys){

    //! Make wenoLevels_ an empty set
    wenoLevels_.clear();

    //! Insert selected keys
    for (auto const& key : keys){
        if (reconstMethods_.count(key) != 0){
            wenoLevels_.insert(key);
        } else {
            cout << "Weno reconstruction level" << key << " has not been defined. " << endl;
            break;
        }
    }
 
}

/**
 * Evaluation of given point with selected weno reconstruction method.
 */
double multiLevelReconstruction::EvaluateMLWENO(const MeshInfo& mi, vertex point, indice localCell){

    double work = 0.0;

    //! Convert local cell indice to global cell indice.
    indice globalCell = mi.MPIlocalCellStart + localCell;

    //! Get calculated nonlinear weights linked to this cell.
    unordered_map<std::string, unordered_map<int,double>> nlw = nonLinearWgts_[FlatIndic(mi,globalCell)];

    //! Evaluate in the multi level weno fashion.
    for (auto const& level : wenoLevels_){
        if (nlw[level].empty() == 0){
            for (auto & wgts: nlw[level]){
                indice owner = localCell + Bend(reconstLevels_[level]->GetSizeX(),wgts.first);
                work += wgts.second * reconstLevels_[level]->Evaluate(mi, owner, point); 
            }
        }
    }

    return work;
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
    cout << "There are " <<reconstLevels_.size()<< " levels pre created." << endl;

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

void multiLevelReconstruction::PrintNonLinearWgts(const MeshInfo& mi){

    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){

        indice add {i,j};
        indice start = mi.MPIlocalCellStart + add;
        unordered_map<std::string, unordered_map<int,double>> nlw = nonLinearWgts_[FlatIndic(mi,start)];

        cout << "Reconstruction at cell ( " << start[0] << ", " << start[1] << ")" << endl; 
        for (auto const& level : wenoLevels_) {
            int sizeX = reconstLevels_[level]->GetSizeX();

            if (nlw[level].empty() == 0) {
                for (auto & in:nlw[level]){
                    indice m = Bend(sizeX,in.first);
                    cout << "At Level " << level << " reconstruction at ( " << m[0] << ", "
                         << m[1] << ") " << " with wgt " << in.second << endl;
                }
            }
        }
    }}
}

void multiLevelReconstruction::Clear(){

    for (auto& it : reconstLevels_){
        delete it.second;
    }

    reconstLevels_.clear();
    reconstMethods_.clear();

    interiorCells_.clear();
    boundaryCells_.clear();

    wenoLevels_.clear();
    interiorLevels_.clear();
    boundaryLevels_.clear();

    nonLinearWgts_.clear();

    etaBias_.clear();
}
