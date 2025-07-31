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

int tensorstencilpoly::setCoef(const MeshInfo& mi,
                               const int& gstartx,   const int& gstarty){

    // Extract tensor product stencil
    vector<vector<vertex>> cornerSet;
    vector<vertex> refcell;

    for (int j=0; j<sizey; j++){
        for (int i=0; i<sizex; i++){
            indice global {i+gstartx, j+gstarty};
            cornerSet.push_back(extractCorners(mi, global));
        }
    }

    // Setup linear system for computing basis polynomials 
    lapack_int n    = sizex*sizey;
    lapack_int nrhs = n;
    lapack_int lda  = n;
    lapack_int ldb  = nrhs;

    coef = new double [n*nrhs] ();

    double * a = new double [n*n] ();
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

    int err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, p, coef, ldb);

    if (err){
        printf("ERROR: Weno Basis Coefficient for order %d, %d. Error type %d \n",
                 sizex,sizey,err);
    }

    delete [] a;
    delete [] p;

    return 1;
}
