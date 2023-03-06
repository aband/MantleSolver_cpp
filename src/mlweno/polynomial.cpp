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

    //! Initialize smoothness Indicator each time it computes
    smoothnessIndic_ = 0.0;

/*
 *    if (Xi_.empty()) {CreateXi_();}
 *
 *    for (int alpha1 = 0; alpha1<maxR_; alpha1++){
 *    for (int alpha2 = 0; alpha2<maxR_-alpha1; alpha2++){
 *        smoothnessIndic_ += Xi_[alpha1]*Xi_[alpha2] -
 *                            (1/((2*alpha1+1)*pow(4,alpha1)))*
 *                            (1/((2*alpha2+1)*pow(4,alpha2)));
 *
 *        double coefindx = 0.0;
 *
 *        if (stencilPolyn_.getI() == maxR_){
 *            coefindx = alpha1 + alpha2*maxR_;
 *        } else {
 *            coefindx = alpha2 + alpha1*maxR_;
 *        }
 *
 *        smoothnessIndic_ *= pow(collapsePolyn_->getCoef(coefindx),2);
 *
 *    }}
 *
 */

    //! Subtracting zero derivative case
    //smoothnessIndic_ -= (Xi_[0]*Xi_[0] - 
    //                    (1/((2*0+1)*pow(4,0)))*
    //                    (1/((2*0+1)*pow(4,0)))) * pow(collapsePolyn_->getCoef(0),2);

    //! Special treatment for tensor product polynomial
    //! No need to create auxiliary Xi term this time

    if (stencilPolyn_.getI()*stencilPolyn_.getJ() != 1){

        for (int all = 1; all<stencilPolyn_.getJ()*stencilPolyn_.getI(); all++){
            int l = all/stencilPolyn_.getI();
            int m = all%stencilPolyn_.getI();

            for (int r=l; r<stencilPolyn_.getJ(); r++){
            for (int s=m; s<stencilPolyn_.getI(); s++){
                smoothnessIndic_ += pow(factorial(r,r-l),2)/(2*(r-l)+1)/pow(4,r-l) * 
                                    pow(factorial(s,s-m),2)/(2*(s-m)+1)/pow(4,s-m); 

                smoothnessIndic_ *= pow(collapsePolyn_->getCoef(FlatIndic(stencilPolyn_.getI(),s,r)),2);
            }}
        }

    }

}

//! Evaluation of auxiliary variable Xi in evaluation of smooth indicator for plain polynomials
void stencilPolynomial::CreateXi_(){
   
    if (stencilPolyn_.getI() > stencilPolyn_.getJ()){
        maxR_ = stencilPolyn_.getI();
        minR_ = stencilPolyn_.getJ();
    } else {
        maxR_ = stencilPolyn_.getJ();
        minR_ = stencilPolyn_.getI();
    }

    Xi_.resize(maxR_,0.0);

    for (int i=0; i<Xi_.size(); i++){
        for (int k=0; k<i+1; k++){
            Xi_[i]  += pow(coef_,2*k)*pow(factorial(i,i-k),2)/(2*(i-k)+1)/pow(4,i-k);
        }
    }

}

double stencilPolynomial::GetSmoothIndic(const MeshInfo& mi, const stencil<indice>& stencilIndice){
    EvalSmoothIndic_(mi, stencilIndice);

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
