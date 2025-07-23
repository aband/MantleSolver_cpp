#include "efficient_stenpoly.h"

stenpoly::stenpoly(const int& sizex,
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

int stenpoly::setCoef(const vector<vector<vertex>>& cornerSet,
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

    Tensor<int> porder = Tensor<int>(2);
    porder.setSize({xdegree, ydegree});

    for (int cell=0; cell<n; cell++){

        vector<vertex> work = cornerSet.at(cell);
 
        for (int ypow = 0; ypow < ydegree; ypow++){
            for (int xpow =0; xpow < xdegree; xpow++){
                int r = porder.getIndex({xpow, ypow});
                a[cell*n+r] = NumIntegralFace(work, {xpow, ypow}, center, scale, basePoly);
            }
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
