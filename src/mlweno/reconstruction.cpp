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

void reconstruction::CreateStencilPolynomials(const indice& start,              const vertex& center,
                                              const vector<indice>& targetCell, const MeshInfo& mi) {

    stencilPolyn_.resize(stencilIndice_.size()); 

    linWgts_.resize(stencilIndice_.size());
    nonLinWgts_.resize(stencilIndice_.size());

    for (int s=0; s<stencilIndice_.size(); s++){
        stencilPolyn_[s] = new stencilPolynomial(start, center, targetCell); 
        stencilPolyn_[s]->SetUpScale(mi,stencilIndice_[s]);
        //stencilPolyn_[s]->SetUpScale(mi,targetCell);
        stencilPolyn_[s]->SetStencilPolynomials(mi,stencilIndice_[s]);
        linWgts_[s] = 1.0;
        nonLinWgts_[s] = 1.0;
    }

}

void reconstruction::AddStencilPolynomials(const indice& start,              const vertex& center,
                                           const vector<indice>& targetCell, const MeshInfo& mi) {
    int currentSize = stencilPolyn_.size();

    assert(currentSize < stencilIndice_.size());

    for (int i=currentSize; i<stencilIndice_.size(); i++){
        stencilPolyn_.push_back(new stencilPolynomial(start, center, targetCell));
        stencilPolyn_[i]->SetUpScale(mi, stencilIndice_[i]);
        stencilPolyn_[i]->SetStencilPolynomials(mi, stencilIndice_[i]);
        linWgts_.push_back(1.0);
        nonLinWgts_.push_back(1.0);
    }

}

void reconstruction::PrintStencils() const {

    cout << "Stencil details ... " << endl;

    for (int s=0; s<stencilIndice_.size(); s++){
        stencil <indice> siNow = stencilIndice_[s];
        for (int j=0; j<siNow.getJ(); j++){
        for (int i=0; i<siNow.getI(); i++){
            indice iNow = siNow(i,j);
            cout << std::setw(2) << "(" << iNow[0] << "," << iNow[1] << ") " ;
        } cout << endl;}
        cout << endl;
    }

    cout << "Stencil polynomials ..." << endl;

    for (int s=0; s<stencilPolyn_.size(); s++){
        cout << endl;
        stencilPolyn_[s]->printCoef();
    }

}

void reconstruction::ComputeSmoothnessIndicatorPolyn_(const MeshInfo& mi) {

    smoothnessIndicPolyn_.resize(stencilIndice_.size());

    for (int s=0; s<stencilIndice_.size(); s++) {
        smoothnessIndicPolyn_[s] = stencilPolyn_[s]->GetSmoothIndic(mi, stencilIndice_[s]);
    }
}

void reconstruction::ComputeNonLinWgts_(const MeshInfo& mi){

    ComputeSmoothnessIndicatorPolyn_(mi);

    if (etaBias_.empty()) {etaBias_.resize(stencilPolyn_.size());std::fill(etaBias_.begin(),etaBias_.end(),0);};
    if (scale_ == -1) {scale_ = stencilPolyn_[0]->GetScale();};

    // Compute non linear WENO weights 
    nonLinWgts_.resize(stencilPolyn_.size());
    double sum = 0.0;
    for (int i=0; i<stencilPolyn_.size(); i++){
        //cout << smoothnessIndicPolyn_[i] << endl;
        nonLinWgts_[i] = linWgts_[i] / pow(smoothnessIndicPolyn_[i] + eps0_*scale_*scale_, 
                                           max(stencilPolyn_[i]->GetOrderX(),stencilPolyn_[i]->GetOrderY())) * 
                                       pow(eps0_*scale_ / (smoothnessIndicPolyn_[i] + eps0_*scale_), etaBias_[i]);

        sum += nonLinWgts_[i];
    }

    std::transform(nonLinWgts_.begin(), nonLinWgts_.end(), nonLinWgts_.begin(), [sum](double x){return x/sum;});

}

double reconstruction::Eval(double x, double y) const{
    assert(nonLinWgts_.empty() == 0);

    double work = 0.0;
    // Already calculated Nonlinear weights
    for (int i=0; i<stencilPolyn_.size(); i++){
        work+=nonLinWgts_[i]*stencilPolyn_[i]->eval(x,y); 
    }

    return work;
}

void reconstruction::PrintSmoothnessIndic() const{

    for (int s=0; s<smoothnessIndicPolyn_.size(); s++){
        cout << smoothnessIndicPolyn_[s] << endl;
    }

}

void reconstruction::PrintNonLinWgts() const{
    for (int s=0; s<nonLinWgts_.size(); s++){
        cout << nonLinWgts_[s] << endl;
    }
}

void reconstruction::Clear() {

    shift_.clear(); 

    stencilSizeMax_[0] = 0;
    stencilSizeMax_[1] = 0;

    stencilNum_ = 0;

    stencilIndice_.clear();

    stencilPolyn_.clear();

    smoothnessIndicPolyn_.clear();

    linWgts_.clear();
    nonLinWgts_.clear();
}

// ======================================================================
singleLevelReconstruction::singleLevelReconstruction(int stencilSizeX, int stencilSizeY){
    stencilSizeX_ = stencilSizeX;
    stencilSizeY_ = stencilSizeY;
 
    // Create stencil of indice used in creating stencil polynomials
    stencilIndice_.SetStencil(stencilSizeX_, stencilSizeY_); 

    for (int j=0; j<stencilSizeY_; j++){
    for (int i=0; i<stencilSizeX_; i++){
        stencilIndice_(i,j) = {i,j};
    }}
}

void singleLevelReconstruction::IdentifyInteriorCell_(const MeshInfo& mi){
    // Should be called each time add a new level to reconstruction
    // For better countability, all stencils
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

vertex singleLevelReconstruction::ComputeStencilCenter_(const MeshInfo& mi, int flat){
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

void singleLevelReconstruction::UpdateSmoothnessIndic_(const MeshInfo& mi){
    for (auto& ind:interior_){
        smoothnessIndic_[ind] = singleLevel_[ind]->GetSmoothIndic(mi,stencilIndice_);
    }
}

const double singleLevelReconstruction::CalculateSmoothnessIndic(const MeshInfo& mi, indice owner) {
    return singleLevel_[FlatIndic(mi, owner)]->GetSmoothIndic(mi,stencilIndice_);
}

void singleLevelReconstruction::ComputeStencilPolyn_(const MeshInfo& mi){
    for (auto& flat: interior_){
        singleLevel_[flat] = new stencilPolynomial(Bend(mi,flat), ComputeStencilCenter_(mi,flat));
        singleLevel_[flat]->SetUpScale(mi,stencilIndice_);
        singleLevel_[flat]->SetStencilPolynomials(mi,stencilIndice_);
    }
}

void singleLevelReconstruction::CreateStencilPolynomials(const MeshInfo& mi){
    IdentifyInteriorCell_(mi);
    ComputeStencilPolyn_(mi);
}

void singleLevelReconstruction::CheckStencilPolynomials(const MeshInfo& mi, indice start){
    singleLevel_[FlatIndic(mi,start)]->printCoef();
}

void singleLevelReconstruction::PrintSmoothnessIndicator(const MeshInfo& mi){

    UpdateSmoothnessIndic_(mi);

    for (auto& ind:interior_){
        cout << smoothnessIndic_[ind] << endl;
    }

}

// ==========================================================================================
// Class of multi level reconstructions, managing information related to smoothness indicator
// and nonlinear weights between reconstruction levels.
// ==========================================================================================
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

void multiLevelReconstruction::UpdateTwoStageNonLinearWgts_(const MeshInfo& mi, indice start){

    UpdateNonLinearWgts_(mi,start);




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
                               pow(sm + scale*scale*eps0_ , max(sizeX, sizeY)) * 
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

// For special treatment on boundary
void multiLevelReconstruction::AddSpecialLinearWgts_(const MeshInfo& mi, const unordered_set<indice>& start){


}

void multiLevelReconstruction::AddSpecialReconstMethods_(const MeshInfo& mi, const unordered_set<indice>& start, vector<vector<indice>> rm){


}

void multiLevelReconstruction::AddReconstMethod(const MeshInfo& mi, unordered_set<indice> start, vector< vector<indice> > rm){
    AddSpecialReconstMethods_(mi,start,rm);
    AddSpecialLinearWgts_(mi,start);
}

void multiLevelReconstruction::GetInfo(){

    cout << "There are " <<allLevels_.size()<< " levels." << endl;

    for (int l=0; l<allLevels_.size(); l++){
        allLevels_[l]->CheckStencils(); 
    }
}

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
