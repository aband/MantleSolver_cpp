#include "tensorstencilpoly.h"

tensorstencilpoly::tensorstencilpoly(const int& inorder){

    order = inorder;

    sizex = inorder+1;
    sizey = inorder+1;
}

tensorstencilpoly::tensorstencilpoly(const int& insizex,
                                     const int& insizey,
                                     const int& inorder){
    order = inorder;
    sizex = insizex;
    sizey = insizey;
}

int tensorstencilpoly::setCoef(const vector<vector<vertex>>& cornerSet,
                               const vertex& center, const double& h){

    // Setup linear system for computing basis polynomials 
    lapack_int n    = sizex*sizey;
    lapack_int nrhs = n;
    lapack_int lda  = n;
    lapack_int ldb  = nrhs;

    double * a = new double [n*n] ();
    coef = new double [n*nrhs] ();
    lapack_int * p = new int [n] ();

    for (int cell=0; cell<n; cell++){

        vector<vertex> work = cornerSet.at(cell);

        for (int j=0; j<sizey; j++){
            for (int i=0; i<sizex; i++){
                a[cell*n + j*sizex+i] = NumIntegralFace(work, {i,j}, center, h, basePoly);
            }
        }
    }

    fill(coef,coef+n*nrhs,0);
    for (int i=0; i<nrhs; i++) {coef[i*n+i]=a[n*i];}

    int err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, p, b, ldb);

    if (err){
        printf("ERROR: Weno Basis Coefficient for order %d, %d. Error type %d \n",
                 xdegree,ydegree,err);
    }

    delete [] a;
    delete [] p;

    return 1;
}
