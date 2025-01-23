#include "stencilpolynomial.h"

stencilpolynomial::stencilpolynomial(const int& sizex,
                                     const int& sizey){
    size[0] = sizex;
    size[1] = sizey;

    // For our tensor product stencil polynomial
    // Stencil size determines the stencilpolynomial order
    tensorpoly.setSize(size);

    for (int i=0; i<tensorpoly.getSize(); i++){
        tensorpoly(i) = polynomial(sizex,sizey);
    }

}

int stencilpolynomial::setCoef(const vector<vector<vertex>>& cornerSet,
                               const vertex& center, const double& scale){
   
    // For clarification, degree x and y are equal to stencil width x and y
    int xdegree = size[0];
    int ydegree = size[1];
    // Setup linear system for computing basis polynomials 
    lapack_int n    = xdegree*ydegree;
    lapack_int nrhs = n;
    lapack_int lda  = n;
    lapack_int ldb  = nrhs;

    double * a = new double [n*n] ();
    double * b = new double [n*nrhs] ();
    lapack_int * p = new int [n] ();

    for (int cell=0; cell<n; cell++){

        vector<vertex> work = cornerSet.at(cell);
  
        for (int r=0; r<n; r++){
            int xpow = r%xdegree;
            int ypow = r/ydegree;
            a[cell*n+r] = NumIntegralFace(work, {xpow, ypow}, center, scale, basePoly);
        }
    }

    fill(b,b+n*nrhs,0);
    for (int i=0; i<nrhs; i++) {b[i*n+i]=a[n*i];}

    int err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, p, b, ldb);

    if (err){
        printf("ERROR: Weno Basis Coefficient for order %d, %d. Error type %d \n",
                 xdegree,ydegree,err);
    }

    for (int p=0; p<n; p++){
        for (int r=0; r<n; r++){
            tensorpoly(p).setCoef(r,b[r*n+p]);
        }
    }

    delete [] a;
    delete [] b;
    delete [] p;

    return 1;
}

int stencilpolynomial::printCoef(){

    cout << "In this " << tensorpoly.getSize(0) << ", " << tensorpoly.getSize(1) << " stencil." << endl;
    for (int i=0; i<tensorpoly.getSize(); i++){
        tensorpoly(i).printCoef();
        cout << endl;
    }

    cout << endl;

    return 1;
}

double stencilpolynomial::eval(const Tensor<double>& sol, const vertex& point, 
                               const vertex& center, const double& h){

    assert(sol.getSize() == tensorpoly.getSize());

    vertex trans = (point-center)/h;
    double work = 0.0;

    for(int c = 0; c<sol.getSize(); c++){
        work += sol(c) * tensorpoly(c).eval(trans);
    }

    return work;
}

double stencilpolynomial::cellsigma(const double& area,
                                    const vector<vertex>& corners, 
                                    const vertex& center, 
                                    const double& scale,
                                    const int& index){

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>& gpf   = GaussPointsFace;

    double work;

    vector<vertex> tmp = corners;
    /*
     *Transform original corner coordinates with given
     *parameter h and center point. If no transform, pass
     *in h=1.0 and center point as (0.0,0.0).
     */
    for (auto & p : tmp){
        p -= center;
        p = p/h; 
    }

    for (int i=0; i<gpf.size(); i++){
        valarray<double> mapped = GaussMapPointsFace(gpf[i],tmp);
        double jac = abs(GaussJacobian(gpf[i],tmp));
        double gw = gwf[i];
 
        Tensor
    }

    return work;
}

int stencilpolynomial::preparesigma(const vector<double>& area,
                                    const vector<vector<vertex>>& cornerSet,
                                    const vertex& center, const double& scale){

    int total = size[0] * size[1];

    sigma.setSize({total, total});

    // For aach cell
    for (int c=0; c<area.size(); c++){
        vector<vertex> corners = cornerSet.at(c);  // Extract four corners
        for (auto & p: corners) {p-=center; p/=h;} // Transform points locally 
        vertex mapped = GaussMapPointsFace(gpf[i],corners); // Map to gauss points
        double jac = abs(GaussJacobian(gpf[i], tmp));
        double gw  = gwf[i];
        
    }

    for (int j=0; j<total; j++){
        Tensor<double> der1 = Tensor<double>(2);
        der1.setSize(size);
        tensorpoly(j).evalDer(size[0]-1, size[1]-1, );
        for (int i=0; i<total; i++){  
            Tensor<double> der2 = Tensor<double>(2);
            der2.setSize(size); 
        }
    }

    return 1;
}
