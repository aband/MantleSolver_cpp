#include "polynomial.h"

using namespace MLWENO;

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

// ================================================================================
