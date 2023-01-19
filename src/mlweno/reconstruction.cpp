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
}

void singleLevelReconstruction::IdentifyInteriorCell_(const MeshInfo& mi){
    // Should be called each time add a new level to reconstruction
    // For better countability, all stencils 
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){ 
       if (j+mi.MPIlocalCellStart[1]+stencilSizeY_-1<mi.MPIglobalCellSize[1] &&
           j+mi.MPIlocalCellStart[0]+stencilSizeX_-1<mi.MPIglobalCellSize[0] ){
               interior_.insert(FlatIndic_(mi,i,j));
           }
    }}
}

void singleLevelReconstruction::ComputeStencilPolyn_(const indice& start, const vertex& center, 
                                                     const vector<indice>& targetCell){
    for (auto& flat: interior_){
        singleLevel_[flat] = new stencilPolynomial(start, center, targetCell);
    }

}

void singleLevelReconstruction::CreateStencilPolynomials(const MeshInfo& mi){
    IdentifyInteriorCell_(mi);
}

void multiLevelReconstruction::AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY){
    allLevels_.push_back(new singleLevelReconstruction(stencilSizeX, stencilSizeY)); 
}

void multiLevelReconstruction::Clear(){
    for (int i=0; i<allLevels_.size(); i++){
        delete allLevels_[i];
    }
    allLevels_.clear();
}
