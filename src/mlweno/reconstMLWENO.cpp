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
/**!
 * Select weno levels from mlweno preparation class.
 * Pointers will be stored in an unordered_map.
 */
void multiLevelReconstruction::SelectWenoReconstLevel(const unordered_set<std::string>& keys,
                                                      const MLWENOPrepare& mlpPtr){

    // Insert selected keys
    // Steal single level reconstruction pointer from MLWENOPrepare class
    for (auto const& key: keys){
        if (mlpPtr.allLevels_.count(key) == 0){
            cout << "This given level " << key << " has not been defined in mlp ..." << endl;
            break;
        } else if (Levels_.count(key) == 1){
            cout << "This given level " << key << " has already been included ..." << endl;
            break;
        } else {
            Levels_.insert(std::make_pair(key, mlpPtr.allLevels_.at(key)));
        }
    }

}

/**!
 * Modify reconst method with given key.
 * Only modify existing reconst method.
 * Please use AddLevel to add reconstruction levels.
 */
void multiLevelReconstruction::ModifyReconstMethod(const std::string& key,
                                                   const vector<indice>& newReconstMethod){
    if(Methods_.count(key) !=0){
        Methods_.erase(key);
    }

    Methods_.insert(std::make_pair(key, newReconstMethod));

}

void multiLevelReconstruction::SetUpLinearWgts(const std::string& key, 
                                               const vector<double>& linwgts){

    if (LinWgts_.count(key) !=0){
        LinWgts_.erase(key);
    }

    LinWgts_.insert(std::make_pair(key,linwgts));
}

void multiLevelReconstruction::UpdateNonLinearWgts(const MeshInfo& mi, 
                                                   const std::string& weightType,
                                                   bool (*assignML)(const indice& globalCell,
                                                                    const MeshInfo& mi)){
    // Update non linear weights for all cells in the target domain
    if (nonLinearWgts_.empty()){
        // Initialize non linear weights with assigned domain.
        // GlobalCells will be selected by assignML function.
        // All local portion of the mesh will be looped through.
        for (int j=mi.MPIlocalCellStart[1] - mi.cellGhostLayerSize; 
                 j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] + mi.cellGhostLayerSize; j++){
        for (int i=mi.MPIlocalCellStart[0] - mi.cellGhostLayerSize; 
                 i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] + mi.cellGhostLayerSize; i++){
            indice globalCell {i,j};
            if (assignML(globalCell,mi)){
                UpdateNonLinearWgtsCell_(mi, FlatIndic(mi, globalCell), weightType);
            }
        }}

    } else {
        // Update non linear weights.
        for (auto const& nlw: nonLinearWgts_){
            UpdateNonLinearWgtsCell_(mi, nlw.first, weightType);
        }
    }

}

// New standard computation of non linear weights
// described in multi-level paper.
// No need of parameter of one-stage or two-stage
/*
void multiLevelReconstruction::UpdateNonLinearWgts(const MeshInfo& mi, 
                                                   bool (*assignML)(const indice& globalCell, 
                                                                    const MeshInfo& mi)){
    // Update non linear weights for all cells in the target domain
        


}
*/

inline int powerShift(const int& r){

    switch (r) {
        case -1:
            return 0;
        case 1:
            return 1; 
        case 2:
            return 3;
        default: 
            return 4;
    }
}

void multiLevelReconstruction::UpdateNonLinearWgtsCell_(const MeshInfo& mi,
                                                        const int& globalCell,
                                                        const std::string& weightType){
    // Incorportated both one stage and two stage weighting scheme.
    // Power shift set default to be 0,0,0 for one stage.
    unordered_map<std::string, unordered_map<int,double>> nlw;

    double sum = 0.0;

    for (auto const& level : Levels_){
        const int sizeX = level.second->GetSizeX();
        const int sizeY = level.second->GetSizeY();

        for (auto const& rm : Methods_.at(level.first)){
            indice owner = Bend(mi, globalCell) + rm; 
            if (level.second->CheckExist(mi,owner)){
                double scale = level.second->GetScale(FlatIndic(mi, owner));
                double sm = level.second->GetSmoothnessIndic(mi, owner);
                int order = max(sizeX, sizeY);
                int r = (weightType=="one_stage") ? order : -1; 

                double value = 1.0/pow(sm + scale*scale*eps0_, order + powerShift(r));
                nlw[level.first].insert(std::make_pair(FlatIndic(sizeX, rm),value));
                sum += value;
            }
        }
    }

    for (auto const& level: Levels_){
        if (nlw[level.first].empty() == 0){
            for (auto & in : nlw.at(level.first)){
                in.second = in.second/sum;
            }
        }
    }

    nonLinearWgts_.erase(globalCell);
    nonLinearWgts_.insert(std::make_pair(globalCell, nlw));
}

// New nonlinear weights defined in mlweno paper
// No distinguish of one-stage and two-stage nonlinear weighting
// The new weighting method strategicly equivalent to two-stage weighting.
void multiLevelReconstruction::UpdateNonLinearWgtsCell_(const MeshInfo& mi,
                                                        const int& globalCell){

    unordered_map<std::string, unordered_map<int, double>> nlw;
  
    double sum = 0.0;

    for (auto const& level : Levels_){
        const int sizeX = level.second->GetSizeX();
        const int sizeY = level.second->GetSizeY();

        // In tensor product polynomial, 
        // sizeX == sizeY always stands
        // ==========================================================================
        int rl = max(sizeX, sizeY);
        int nl = 1;

        double s = 1;

        if (rl == 1){
            nl = 1;
        } else if (rl == 2){
            nl = 3;
        } else {
            nl = 4;
        }

        //for (auto const& rm: Methods_.at(level.first)){
        for (int k=0; k<Methods_.at(level.first).size(); k++){
            indice rm  = Methods_.at(level.first).at(k);
            indice owner = Bend(mi, globalCell) + rm;
            // Extract linear weight
            double omega_l = LinWgts_.at(level.first).at(k);

            if (level.second->CheckExist(mi, owner)){
                // Scale factor associated with each cells
                double h0 = level.second->GetScale(FlatIndic(mi, owner));
                // Smoothness indicator associated with each cell
                double dm = level.second->GetSmoothnessIndic(mi, owner);

                double omega_hat = omega_l/pow(dm+eps0_*h0*h0,s*rl+nl); 

                nlw[level.first].insert(
                    std::make_pair(FlatIndic(sizeX, rm),omega_hat));
                sum += omega_hat;
            }
        }
    }

    for (auto const& level: Levels_) {
        if (nlw[level.first].empty() == 0){
            for (auto & in : nlw.at(level.first)){
                in.second = in.second/sum;
            }
        }
    }

    nonLinearWgts_.erase(globalCell);
    nonLinearWgts_.insert(std::make_pair(globalCell, nlw));
}

double multiLevelReconstruction::EvaluateMLWENO (const MeshInfo& mi,
                                                 const vertex& point,
                                                 const indice& globalCell) const {
    double work = 0.0;

    unordered_map<std::string, unordered_map<int, double>> nlw = nonLinearWgts_.at(FlatIndic(mi,globalCell));

    for (auto const& level : nlw){
        if (nlw[level.first].empty() == 0){
            for (auto const& wgts : level.second){
                indice owner = globalCell + Bend(Levels_.at(level.first)->GetSizeX(), wgts.first);

                work += wgts.second * Levels_.at(level.first)->Evaluate(mi, owner, point);
            }
        }
    }

    return work; 
}

void multiLevelReconstruction::PrintNonLinearWgts(const MeshInfo& mi){

    for (int j=mi.MPIlocalCellStart[1] - mi.cellGhostLayerSize; 
             j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] + mi.cellGhostLayerSize; j++){
    for (int i=mi.MPIlocalCellStart[0] - mi.cellGhostLayerSize; 
             i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] + mi.cellGhostLayerSize; i++){

        if (i>-1 && i<mi.MPIglobalCellSize[0] &&
            j>-1 && j<mi.MPIglobalCellSize[1]){
            indice start {i,j};
            unordered_map<std::string, unordered_map<int,double>> nlw = nonLinearWgts_[FlatIndic(mi,start)];

            cout << "Reconstruction at cell ( " << start[0] << ", " << start[1] << ")" << endl; 
            for (auto const& level : Levels_) {
                int sizeX = level.second->GetSizeX();
    
                if (nlw[level.first].empty() == 0) {
                    for (auto & in:nlw[level.first]){
                        indice m = Bend(sizeX,in.first);
                        cout << "At Level " << level.first << " reconstruction at ( " << m[0] << ", "
                             << m[1] << ") " << " with wgt " << in.second << endl;
                    }
                }
            }
        }
    }}
}
