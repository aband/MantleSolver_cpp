#include "polynomial.h"

using namespace MLWENO;

double basePoly(vertex& point, const vector<int>& param){
    return pow(point[0],param[0])*pow(point[1],param[1]);
}

// Constructors and destructors
basisPolynomial::basisPolynomial(const int maxDegree[2]){
    maxDegree_[0] = maxDegree[0];
    maxDegree_[1] = maxDegree[1];

    //if (coef_) {delete [] coef_;}
    coef_ = new double [maxDegree[0]*maxDegree[1]];
}

basisPolynomial::basisPolynomial(const int maxDegree[2], double* coef){
    maxDegree_[0] = maxDegree[0];
    maxDegree_[1] = maxDegree[1];

    //if (coef_) {delete [] coef_;}
    int coefSize = maxDegree[0]*maxDegree[1];
    coef_ = new double [coefSize]();
    for (int i=0; i<coefSize; i++){coef_[i] = coef[i];}

}

basisPolynomial::~basisPolynomial(){
    delete [] coef_;
}

void basisPolynomial::setMaxDegree(const int maxDegree[2]){
    maxDegree_[0] = maxDegree[0];
    maxDegree_[1] = maxDegree[1];
}

void basisPolynomial::setCoef(double* coef){

    if (maxDegree_[0] == -1 || maxDegree_[1] == -1) {
        cout <<" Max Degrees weren't assigned ! " << endl;
    }

    if (coef_) {delete [] coef_;}
    int coefSize = maxDegree_[0]*maxDegree_[1];
    coef_ = new double [coefSize]();
    for (int i=0; i<coefSize; i++){coef_[i] = coef[i];}
}

double basisPolynomial::eval(double x, double y) const {
    // Evaluation of the 2D basis polynomial with the given point
    // Using Horner's method
    double ycoef[maxDegree_[1]];

    int start = 0;

    for (int r=0; r<maxDegree_[1]; r++ ){
        ycoef[r] = polyEval(x,&coef_[start],maxDegree_[0]-1);
        start += maxDegree_[0];
    }

    return polyEval(y,ycoef,maxDegree_[1]-1);
}

double* basisPolynomial::getCoef() const{
    double * coef = new double [maxDegree_[0]*maxDegree_[1]] ();
    for (int i=0; i<maxDegree_[0]*maxDegree_[1]; i++){
        coef[i] = coef_[i];
    }
    return coef;
} 

double basisPolynomial::getCoef(int i) const{
    assert(coef_ != nullptr); 

    return coef_[i];  

}

void basisPolynomial::printCoef() const {
    for (int i=0; i<maxDegree_[0]*maxDegree_[1]; i++){
        cout << std::setprecision(5)<< coef_[i] << "  " ;
    }cout << endl;
}

// ================================================================================
stencilPolynomial::stencilPolynomial(const indice& start, const vertex& center){

    start_[0] = start[0];
    start_[1] = start[1];

    center_[0] = center[0];
    center_[1] = center[1];

}

stencilPolynomial::stencilPolynomial(const indice& start, const vertex& center, 
                                     const vector<indice>& targetCell){


    start_[0] = start[0];
    start_[1] = start[1];

    center_[0] = center[0];
    center_[1] = center[1];

    targetCell_.resize(targetCell.size());

    for (int i=0; i<targetCell_.size(); i++){
        targetCell_.at(i) = targetCell.at(i);
    }

}

void stencilPolynomial::SetTargetCell_(const vector<indice>& targetCell){
    targetCell_.resize(targetCell.size());

    for (int i=0; i<targetCell_.size(); i++){
        targetCell_.at(i) = targetCell.at(i);
    }
}

void stencilPolynomial::ComputeCellBasedScale_(const MeshInfo& mi){
    scale_ = 0.0;
    for (auto & cell: targetCell_){
        vector<vertex> work;
        for (auto & c: mi.faceCorner){
            int sj = start_[1]+cell[1]+c[1] + mi.vertexGhostLayerSize;
            int si = start_[0]+cell[0]+c[0] + mi.vertexGhostLayerSize;
            work.push_back(mi.lmesh[sj*mi.MPIlocalVertexSizeFull.at(0)+si]);
        }
        scale_ += NumIntegralFace(work, {0,0}, {0.0,0.0}, 1.0, constFunc); 
    }

    scale_ = pow(scale_,0.5); 
}

void stencilPolynomial::ComputeStencilBasedScale_(const MeshInfo& mi, const stencil <indice>& stencilIndice){

    double maxScale_ = 0.0;

    for (int j=0; j<stencilIndice.getSize(); j++){
        scale_ = 0.0; 
        vector<vertex> work;
        indice siNow = stencilIndice(j);
        for (auto & c: mi.faceCorner){
            int sj = start_[1]+siNow[1]+c[1] + mi.vertexGhostLayerSize;
            int si = start_[0]+siNow[0]+c[0] + mi.vertexGhostLayerSize;
            work.push_back(mi.lmesh[sj*mi.MPIlocalVertexSizeFull.at(0)+si]);
        }
        scale_ += NumIntegralFace(work, {0,0}, {0.0,0.0}, 1.0, constFunc); 

        if (scale_ > maxScale_) {maxScale_ = scale_;};

    }

    scale_ = pow(maxScale_,0.5);

}

void stencilPolynomial::SetStencilPolynomials(const MeshInfo& mi, 
                                              const stencil <indice>& stencilIndice){

    stencil<indice> siNow = stencilIndice;

    // Setup linear system for computing basis polynomials 
    lapack_int n    = siNow.getSize();
    lapack_int nrhs = n;
    lapack_int lda  = n;
    lapack_int ldb  = nrhs;

    double * a = new double [n*n] ();
    double * b = new double [n*nrhs] ();
    lapack_int * p = new int [n] ();

    for (int cell = 0; cell<n; cell++){
        // Cell indice  
        indice currentCell = start_ + siNow(cell);

        vector<vertex> work;
        for (auto & c: mi.faceCorner){
            int sj = currentCell[1] + c[1] + mi.vertexGhostLayerSize;
            int si = currentCell[0] + c[0] + mi.vertexGhostLayerSize;

            work.push_back(mi.lmesh[sj*mi.MPIlocalVertexSizeFull.at(0)+si]);
        }

        for (int r = 0; r<n; r++){
            int xpow = r%siNow.getI();
            int ypow = r/siNow.getI();
            a[cell*n + r] = NumIntegralFace(work, {xpow,ypow}, center_, scale_, basePoly);
        }
    }
    fill(b,b+n*nrhs,0);
    for (int i=0; i<nrhs; i++) {b[i*n+i]=a[n*i];}

    int err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, p, b, ldb);

    int maxDegree[2] = {siNow.getI(),siNow.getJ()};

    if (err){
        printf("ERROR: Weno Basis Coefficient for order %d, %d. Error type %d \n",
                 maxDegree[0],maxDegree[1],err);
    }

    stencilPolyn_.SetStencil(siNow.getI(),siNow.getJ());

    for (int p=0; p<n; p++){
        double * tmpcoef = new double [n]();
       
        for (int r=0; r<n; r++){
            tmpcoef[r] = b[r*n+p];
        }

        stencilPolyn_(p) = new basisPolynomial(maxDegree, tmpcoef);

        delete [] tmpcoef;
    }

    delete [] a;
    delete [] b;
    delete [] p;

}

void stencilPolynomial::SetCollapsePolyn(const MeshInfo& mi, const stencil <indice>& stencilIndice) {
    SetCollapsePolyn_(mi,stencilIndice);
}

void stencilPolynomial::SetCollapsePolyn_(const MeshInfo& mi, const stencil <indice>& stencilIndice) {
    stencil<indice> siNow = stencilIndice;
    double * tmpcoef = new double [stencilPolyn_.getSize()]();
    for (int i=0; i<stencilPolyn_.getSize(); i++){
        double * tmp = stencilPolyn_(i)->getCoef();
        indice currentCell = start_ + siNow(i);

        for (int j=0; j<stencilPolyn_.getSize(); j++){
            tmpcoef[j] += tmp[j]*mi.localVals[currentCell[1]][currentCell[0]];
        }
        delete [] tmp;
    }
    int maxDegree[2] = {stencilPolyn_.getI(),stencilPolyn_.getJ()};
    collapsePolyn_ = new basisPolynomial(maxDegree,tmpcoef);
    delete [] tmpcoef;
} 

void stencilPolynomial::EvalSmoothIndic_(const MeshInfo& mi, const stencil <indice>& stencilIndice){

    SetCollapsePolyn_(mi, stencilIndice);

    // Initialize smoothness Indicator
    smoothnessIndic_ = 0.0;

    // ===================================

    for (int i=1; i<collapsePolyn_->getSize(); i++){
        int l = i/collapsePolyn_->getI(); int m = i%collapsePolyn_->getI(); 
        for (int j=i; j<collapsePolyn_->getSize(); j++){
            int r = j/collapsePolyn_->getI(); int s = j%collapsePolyn_->getI();
            smoothnessIndic_ += pow(collapsePolyn_->getCoef(j) *
                                //factorial(r,r-l)*factorial(s,s-m),2)*
                                (factorial(r)/factorial(l)) * (factorial(s)/factorial(m)),2) *
                                pow(2*(double)(r-l)+1.0,-1)*
                                pow(2*(double)(s-m)+1.0,-1);

        }
    }
}

double stencilPolynomial::GetSmoothIndic(const MeshInfo& mi, const stencil<indice>& stencilIndice){
    if (smoothnessIndic_ == -1) {
        EvalSmoothIndic_(mi, stencilIndice);
    }

    return smoothnessIndic_;
}

double stencilPolynomial::eval(double x, double y) const{

    assert(collapsePolyn_ != nullptr);

    return collapsePolyn_->eval(x,y);
}

void stencilPolynomial::printCoef() {
    for (int s =0; s<stencilPolyn_.getSize(); s++){
        stencilPolyn_(s)->printCoef();
    }
}

void stencilPolynomial::printCoef(int s) {
    assert(s < stencilPolyn_.getSize());
    stencilPolyn_(s)->printCoef();
}
