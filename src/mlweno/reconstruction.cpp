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

void singleLevelReconstruction::UpdateDerivSmoothnessIndic(const MeshInfo& mi){
    for (auto& ind:interior_){
        smoothnessIndicDeriv_.at(ind) = singleLevel_.at(ind)->GetDerivSmoothIndic(mi,stencilIndice_);
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
double singleLevelReconstruction::Evaluate(const MeshInfo& mi, indice owner, vertex point){

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
    //! Print added levels and reconstruction methods
    cout << "There are " <<allLevels_.size()<< " levels created." << endl;

    for (auto const& it : allLevels_){
        cout << it.first << " " ; 
        (it.second)->CheckStencils();
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

//! A more advanced boudarylayer separation function
//! Assign not only the (1,1) level to boundary cells.
void multiLevelReconstruction::SeparateBoundaryLayer(const MeshInfo& mi, const int& layerSize, 
                                                     const unordered_set<std::string>& additionalLevels){

    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){
        indice shift {i,j}; 
        indice global = mi.MPIlocalCellStart + shift;
        if ((global[0] < 0 + layerSize || global[0] > mi.MPIglobalCellSize[0] - layerSize -1) ||
            (global[1] < 0 + layerSize || global[1] > mi.MPIglobalCellSize[1] - layerSize -1)){
            boundaryCells_.insert(FlatIndic(mi,global));
        } else {
            interiorCells_.insert(FlatIndic(mi,global));
        }
    }}

    //! Create two different weno reconstruction levels for interior and boundary cells.
    SeparateReconstMethods_(additionalLevels);
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

void multiLevelReconstruction::SeparateReconstMethods_(const unordered_set<std::string>& additionalLevels){
    interiorLevels_ = wenoLevels_;
    boundaryLevels_ = wenoLevels_;
    for (const auto& al : additionalLevels){
        interiorLevels_.erase(al);
    }
    assert(interiorLevels_.empty() == 0);
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
                if (sizeX*sizeY == 1){power = 1; linearWgt = 1e-2;} else {power = 2; linearWgt = 1;};

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

                double value = wgts[FlatIndic(sizeX,rm)]/pow(sm + scale*scale*eps0_ , max(sizeX, sizeY));
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

void multiLevelReconstruction::UpdateNonLinearWgtsAndDerivs_(const MeshInfo& mi, int flatGlobal, 
                                                             const unordered_set<std::string>& levels){

// This function update non linear weights and the corresponding derivatives at
// the same time!!! Only used when full differentiation is required.

    unordered_map<std::string, unordered_map<int, double>> tmp;
    unordered_map<std::string, unordered_map<int, double>> nlw;

    unordered_map<std::string, unordered_map<int, derivative>> tmpdnlw;
    unordered_map<std::string, unordered_map<int, derivative>> dnlw;

    double sum = 0.0;
    derivative sumdnlw;

    // First stage of the calculation of non linear weights
    // and its derivatives.
    for (auto const& level : levels){
        const int sizeX = reconstLevels_[level]->GetSizeX();
        const int sizeY = reconstLevels_[level]->GetSizeY();
        for (auto const& rm : reconstMethods_[level]){
            indice owner = Bend(mi,flatGlobal) + rm;
            if (reconstLevels_[level]->CheckExist(mi,owner)){
                double scale = reconstLevels_[level]->GetScale(FlatIndic(mi, owner));
                // Get smoothness indicator for each level
                double sm = reconstLevels_[level]->GetSmoothnessIndic(mi, owner);
                // Get derivative of smoothness indicator for each level
                derivative dsm = reconstLevels_[level]->GetSmoothnessIndicDeriv(mi, owner);

                int power = 2;
                double linearWgt = 1.0;
                if (sizeX*sizeY == 1){power = 1; linearWgt = 1e-2;} else {power = 2; linearWgt = 1;};

                // Calculate stage 1 scale value
                double value = linearWgt/pow(sm + scale*scale*eps0_ , power);
                // Calculate derivative of stage 1 scale value
                double modify = linearWgt * (-1*power)/ pow(sm + scale*scale*eps0_, power+1);
                unordered_map_arithmetic(dsm, modify, std::multiplies<double>());

                // Insert stage 1 value into temporary storage
                tmp[level].insert(std::pair<int, double>(FlatIndic(sizeX,rm), value));

                tmpdnlw[level].insert(std::pair<int, unordered_map<int,double>>(FlatIndic(sizeX,rm),dsm));

                // Create sum of all smoothness indicator values
                sum += value;

                unordered_map_arithmetic(sumdnlw, dsm, std::plus<double>());
            }
        }
    }

    // Create stage 1 nonlinear weights
    for (auto const& level: levels){
        if (tmp[level].empty() == 0){
            for (auto & in: tmp[level]){
                in.second = in.second/sum;
            }
        }
    }

    // Create derivative of stage 1 nonlinear weights
    for (auto const& level: levels){
        if (tmpdnlw[level].empty() == 0){
            for (auto& in: tmpdnlw[level]){
                unordered_map_arithmetic(in.second,sum,std::divides<double>());
                double modify = -1 /sum/sum * tmp[level].at(in.first);
                unordered_map_arithmetic(sumdnlw,modify,std::multiplies<double>());
                unordered_map_arithmetic(in.second, sumdnlw, std::plus<double>());
            }
        }
    }

    //! Second stage of the calculatino of non linear weights
    sum = 0.0; 
    sumdnlw.clear();

    for (auto const& level : levels){
        const int sizeX = reconstLevels_[level]->GetSizeX();
        const int sizeY = reconstLevels_[level]->GetSizeY();
        const unordered_map<int,double>& wgts = tmp.at(level);
        const unordered_map<int,derivative>& dwgts = tmpdnlw.at(level);

        for (auto const& rm : reconstMethods_[level]){
            indice owner = Bend(mi,flatGlobal) + rm;
            if (reconstLevels_[level]->CheckExist(mi,owner)){
                double scale = reconstLevels_[level]->GetScale(FlatIndic(mi, owner));
                // Get smoothness indicator for each level
                double sm = reconstLevels_[level]->GetSmoothnessIndic(mi, owner);
                // Get derivative of smoothness indicator for each level
                derivative dsm = reconstLevels_[level]->GetSmoothnessIndicDeriv(mi, owner);

                double modify = pow(sm + scale*scale*eps0_ , max(sizeX, sizeY));

                // Update nonlinear weights
                double value = wgts.at(FlatIndic(sizeX,rm))/modify;
                nlw[level].insert(std::pair<int, double>(FlatIndic(sizeX,rm), value));
                sum += value;

                // Update derivatives of nonlinear weights
                derivative dstage2 = dwgts.at(FlatIndic(sizeX,rm));
                unordered_map_arithmetic(dstage2,modify,std::multiplies<double>());
                double modify2 = value * -1*max(sizeX,sizeY)/pow(sm+scale*scale*eps0_,max(sizeX,sizeY)+1);
                unordered_map_arithmetic(dsm,modify2,std::multiplies<double>());
                unordered_map_arithmetic(dstage2,dsm,std::plus<double>());

                dnlw[level].insert(std::pair<int, derivative>(FlatIndic(sizeX,rm),dstage2));
                unordered_map_arithmetic(sumdnlw,dstage2,std::plus<double>());
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

    // Create derivative of two stage non linear weights
    for (auto const& level: levels){
        if (dnlw[level].empty() == 0){
            for (auto& it: dnlw[level]){
                unordered_map_arithmetic(it.second,sum,std::divides<double>());
                double modify = -1 /sum/sum * nlw[level].at(it.first);
                unordered_map_arithmetic(sumdnlw,modify,std::multiplies<double>());
                unordered_map_arithmetic(it.second, sumdnlw, std::plus<double>());
            }
        }
    }

    derivNonLinearWgts_.erase(flatGlobal);
    derivNonLinearWgts_.insert(std::pair<int, unordered_map<std::string, unordered_map<int,derivative>>>(flatGlobal,dnlw));
 
}

void multiLevelReconstruction::UpdateNonLinearWgts(const MeshInfo& mi){
//! Default multi level reconstruction will utilize two stage nonlinear weights
        for (auto const& it: interiorCells_){
            UpdateTwoStageNonLinearWgts_(mi, it, interiorLevels_);
        }

        for (auto const& it: boundaryCells_){
            UpdateTwoStageNonLinearWgts_(mi, it, boundaryLevels_);
        }
}

void multiLevelReconstruction::UpdateNonLinearWgts(const MeshInfo& mi, const int stage){

    if (prepare_ == false){

        for (auto const& singleLevel : wenoLevels_){
            //cout << "Level " << singleLevel << " Smoothness indic updated ... " << endl;
            reconstLevels_[singleLevel]->UpdateSmoothnessIndic(mi);
        }
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
    // Make sure it is a existing key
    //assert(reconstMethods_.count(key) != 0); 

    if (reconstMethods_.count(key) != 0){
        // Corresponding reconstruction method has been defined beforehand.
        // Ereasing the existing key
        reconstMethods_.erase(key);
    }

        // Assign key and reconstruction method pair
        reconstMethods_.insert(std::pair<std::string, vector<indice>>(key, newReconstMethod));
}

/**
 * Modify wenoLevels_.
 * Selecting predefined reconstruction method in the actual computation process.
 */
void multiLevelReconstruction::SelectWenoReconstLevel(const unordered_set<std::string>& keys){

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
 * "Stealing" prepared single level reconstruction pointer from MLWENOPrepare class.
 */
void multiLevelReconstruction::SelectWenoReconstLevel(const unordered_set<std::string>& keys,
                                                      const MLWENOPrepare& mlpPtr){

    // Indicating MLWENOPrepare is used in the computation.
    prepare_ = true;

    // Make wenoLevels_ an empty set
    wenoLevels_.clear();

    // Insert selected keys
    for (auto const& key: keys){
        if (mlpPtr.allLevels_.count(key) == 0){
            cout << "The given level " << key << " has not been defined ...";
            break;
        } else {
            wenoLevels_.insert(key);
            reconstLevels_[key] = mlpPtr.allLevels_.at(key);
        }
    }

    lowestLevel_ = reconstLevels_.begin()->first;
    highestLevel_ = reconstLevels_.rbegin()->first;
}

/**
 * Reconstruction of a given point value with selected weno reconstruction method.
 */
double multiLevelReconstruction::EvaluateMLWENO (const MeshInfo& mi, vertex point, indice globalCell) const {

    double work = 0.0;

    //! Convert local cell indice to global cell indice.
    //indice globalCell = mi.MPIlocalCellStart + localCell;

    //! Get calculated nonlinear weights linked to this cell.
    unordered_map<std::string, unordered_map<int,double>> nlw = nonLinearWgts_.at(FlatIndic(mi,globalCell));

    //! Evaluate in the multi level weno fashion.
    for (auto const& level : wenoLevels_){
        if (nlw[level].empty() == 0){
            for (auto & wgts: nlw[level]){
                indice owner = globalCell + Bend(reconstLevels_.at(level)->GetSizeX(),wgts.first);
                work += wgts.second * reconstLevels_.at(level)->Evaluate(mi, owner, point); 
            }
        }
    }

    return work;
}

/**
 * Derivative of the reconstructionof value with respect to the given point
 * with multi level weno method.
 * Start with pseudo derivative where non linear weights are not differentiated.
 */
unordered_map<int, double> multiLevelReconstruction::EvaluateDerivMLWENO(const MeshInfo& mi,
                                                                         const vertex& point, 
                                                                         const indice& global) const {

    unordered_map<int, double> work;
    
    //! Extract nonliear weights linked to the target cell.
    unordered_map<std::string, unordered_map<int, double>> nlw = nonLinearWgts_.at(FlatIndic(mi,global));

    //! Evaluate in the multi level weno fashion.
    for (auto const& level : wenoLevels_){
        if (nlw[level].empty() == 0){
            for (auto & wgts: nlw[level]){
                indice owner = global + Bend(reconstLevels_.at(level)->GetSizeX(), wgts.first);
                int flatOwner = FlatIndic(mi, owner);

                for (int p = 0; p<reconstLevels_.at(level)->GetSizeX()*
                                  reconstLevels_.at(level)->GetSizeY(); p++){
                    indice ownerShift = owner + Bend(reconstLevels_.at(level)->GetSizeX(),p);
                    int flatOwnerShift = FlatIndic(mi,ownerShift);

                    if (work.count(flatOwnerShift)>0){
                        work[flatOwnerShift] += wgts.second * 
                        reconstLevels_.at(level)->Evaluate(mi,owner,point, p);
                    } else {
                        work.insert(std::pair<int, double>(flatOwnerShift, wgts.second* 
                                    reconstLevels_.at(level)->Evaluate(mi,owner,point,p)));
                    }
                }
            }
        }
    }

    return work;
}

/**
 * Differentiate non linear weight when 
 * evaluate MLWENO derivative.
 */
unordered_map<int, double> multiLevelReconstruction::EvaluateDerivMLWENO(const MeshInfo& mi,
                                                     const vertex& point, const indice& global,
                                                     const int& flag) const{
    assert(flag == 1);

    unordered_map<int, double> work;

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

//! Clear created multi level reconstruction
void multiLevelReconstruction::Clear(){

    if (prepare_ == 0){
        // Clear pointers if MLWENOPrepare is not used
        // Otherwise, these pointers should be destroyed in MLWENOPrepare
        for (auto& it : reconstLevels_){
            delete it.second;
        }
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
