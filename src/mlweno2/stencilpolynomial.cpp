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

double stencilpolynomial::sigmaintegral(const vector<vertex>& corners, 
                     const double& area,
                     const vertex& center, const double& h,
                     const int& index1, const int& index2){

    double work = 0.0;

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    vector<vertex> tmp = corners;          // Extract four corners
    for (auto & p: tmp) {p-=center; p/=h;} // Transform points locally

    int total = size[0] * size[1];  

    Tensor<double> der = Tensor<double>(2);
    der.setSize(size);
    Tensor_zero(der);

    for (int i=0; i<gpf.size(); i++){
        valarray<double> mapped = GaussMapPointsFace(gpf[i],tmp);
        double jac = abs(GaussJacobian(gpf[i],tmp));
        double gw = gwf[i];

        Tensor<double> der1 = Tensor<double>(2);
        der1.setSize(size);
        Tensor<double> der2 = Tensor<double>(2);
        der2.setSize(size);

        tensorpoly(index1).evalDer(size[0]-1,size[1]-1,mapped[0], mapped[1], 1.0, 
                                   der1);
        tensorpoly(index2).evalDer(size[0]-1,size[1]-1,mapped[0], mapped[1], 1.0, 
                                   der2);

        Tensor_multi_add(der1,der2,gw*jac,der);
    }

    // Sum through all order of derivatives
    for (int i=0; i<total; i++){
        work += der(i) * pow(area/h*h,i); 
    }

    return work;
}

int stencilpolynomial::sigma(const vector<vertex>& corners, const double& area,
                             const vertex& center, const double& h){

    int total = size[0] * size[1];
    tensorsigma.setSize({total, total});

    // diagonal first
    for (int d=0; d<total; d++)
    {tensorsigma({d,d}) = sigmaintegral(corners, area, center, h, d, d);}
 
    // symmetrical applied by computing upper triangle only
    for (int j=1; j<total; j++){
        for (int i=0; i<j; i++){
            tensorsigma({i,j}) = sigmaintegral(corners, area, center, h, i, j);
            tensorsigma({j,i}) = tensorsigma({i,j});
        }
    }

    return 1;
}
