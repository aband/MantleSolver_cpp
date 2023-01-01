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

// ================================================================================
stencilPolynomial::stencilPolynomial(const indice& start, const vertex& center, 
                                     const vector<indice>& targetCell, const MeshInfo& mi){

    start_[0] = start[0];
    start_[1] = start[1];

    center_[0] = center[0];
    center_[1] = center[1];

    targetCell_.resize(targetCell.size());

    for (int i=0; i<targetCell_.size(); i++){
        targetCell_.at(i) = targetCell.at(i);
    }

    // Compute scale
    double scale = 0.0;
    for (auto & cell: targetCell){
        vector<vertex> work;
        for (auto & c: mi.faceCorner){
            int sj = start[1]+cell[1]+c[1];
            int si = start[0]+cell[0]+c[0];
            work.push_back(mi.lmesh[sj*mi.MPIlocalSizeFull.at(0)+si]);
        }
    }

}

void stencilPolynomial::SetStencilPolynomials(const MeshInfo& mi, 
                                              const vector<stencil <indice>>& stencilIndice){

    stencilPolyn_.resize(stencilIndice.size()); 

    // Iteration through stencils of indices
    for (int s=0; s<stencilIndice.size(); s++) {

        stencil <indice> siNow = stencilIndice[s];

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
            indice currentCell = start_ + siNow(n);
            for (int r = 0; r<n; r++){
                vector<vertex> work;
                for (auto & c: mi.faceCorner){
                    int sj = currentCell[1] + c[1];
                    int si = currentCell[0] + c[0];

                    work.push_back(mi.lmesh[sj*mi.MPIlocalSizeFull.at(0)+si]);
                }
                int xpow = r%stencilIndice[s].getI();
                int ypow = r/stencilIndice[s].getI();
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

        stencil <basisPolynomial*> singleStencilPolyn(n);

        for (int p=0; p<n; p++){
            double * tmpcoef = new double [n]();
           
            for (int r=0; r<n; r++){
                tmpcoef[r] = b[r*n+p];
            }

            singleStencilPolyn(p) = new basisPolynomial(maxDegree, tmpcoef);

            delete [] tmpcoef;
        }

        delete [] a;
        delete [] b;
        delete [] p;
    }

}
