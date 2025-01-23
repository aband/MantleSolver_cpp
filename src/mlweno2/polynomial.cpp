#include "polynomial.h"

// Evaluate derivative of 1D polynomial up to a given derivative with Horner's method
// Input der starts from 0.
int computeDerivative(const int& der,  const int& degree, 
                      const double& x, const double& scale,
                      double * coef, double * work){

    double xx = x/scale;

    for (int d=0; d<=der; d++){work[d] = 0.0;}

    for (int i=degree-1; i>=0; i--){
        for (int d=der; d>=1; d--){
            work[d] = work[d]*xx + d*work[d-1];
        }
        work[0] = work[0]*xx + coef[i];
    }

    for(int d=1; d<=der; d++) {work[d] /= pow(scale,d);}

    return 1;
}

// Numerical integral function but used specifically for polynomial
double polyNumIntegralFace(const vector<vertex>& corners,
                           const double& h,
                           const vertex& center,
                           polynomial& mypoly){

    // Copy gauss weights and gauss points
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>& gpf = GaussPointsFace;

    assert(corners.size() == 4);

    vector< valarray<double> > tmp = corners;
    /*
     *Transform original corner coordinates with given
     *parameter h and center point. If no transform, pass
     *in h=1.0 and center point as (0.0,0.0).
     */
    for (auto & p : tmp){
        p -= center;
        p = p/h; 
    }

    double work = 0.0;

    for (size_t i=0; i<gpf.size(); i++){
        valarray<double> mapped = GaussMapPointsFace(gpf[i],tmp);
        double jac = abs(GaussJacobian(gpf[i],tmp));
        double gw = gwf[i];
        work += jac*gw*mypoly.eval(mapped[0], mapped[1]); 
    }

    return work;
}

// ===================================================================================

polynomial::polynomial(const int& degreex,
                       const int& degreey){

    degree[0] = degreex;
    degree[1] = degreey;
}

polynomial::~polynomial(){
    delete [] coef;
}

// Only used for testing purpose
int polynomial::resetDegree(const int& degreex, 
                            const int& degreey){

    degree[0] = degreex;
    degree[1] = degreey;

    return 1;
}

int polynomial::setCoef(double* setcoef){

    if (degree[0] == -1 || degree[1] == -1) {
        cout <<" Max Degrees weren't assigned properly ! " << endl;
    }

    if (coef) {delete [] coef;}
    int coefSize = degree[0]*degree[1];
    coef = new double [coefSize]();
    for (int i=0; i<coefSize; i++){coef[i] = setcoef[i];}

    return 1;
}

int polynomial::setCoef(int index, double val){

    if (coef == nullptr){coef = new double [degree[0]*degree[1]]();};

    coef[index] = val;

    return 1;
}

int polynomial::printCoef() const {
    for (int i=0; i<degree[0]*degree[1]; i++){
        cout << std::setprecision(5)<< coef[i] << "  " ;
    }cout << endl;
    return 1;
}

double polynomial::eval(const double& x,
                        const double& y) const{

    double ycoef[degree[1]];

    int start = 0;

    for (int r=0; r<degree[1]; r++ ){
        ycoef[r] = polyEval(x,&coef[start],degree[0]-1);
        start += degree[0];
    }

    return polyEval(y,ycoef,degree[1]-1);
}

int polynomial::evalDer(const int& derx,    const int& dery,
                        const double& x,    const double& y,
                        const double& scale, Tensor<double>& derTensor){

    int degreex = degree[0];
    int degreey = degree[1];

    double * workx = new double [derx + 1];  
    double * worky = new double [dery + 1];

    double ycoef[derx + 1][degreey];

    int start = 0;

    // Horner's method in y
    for (int r=0; r<degree[1]; r++){
        computeDerivative(derx, degreex, x, scale, &coef[start], workx); 
        for (int d=0; d<=derx; d++){ycoef[d][r] = workx[d];}
        start += degreex;
    }

    // Horner's method in x
    for (int dx=0; dx<=derx; dx++){
        computeDerivative(dery, degreey, y, scale, ycoef[dx], worky);
        for (int dy=0; dy<=dery; dy++){
            derTensor({dx,dy}) = worky[dy];
        }
    }

    return 1;
}
