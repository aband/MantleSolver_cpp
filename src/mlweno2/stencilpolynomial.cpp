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

int stencilpolynomial::setCoef(const MeshInfo& mi,
                               const Tensor<indice>& stencilindice,
                               const vertex& center, const double& scale){

    Tensor<indice> siNow = stencilindice;

    // Setup linear system for computing basis polynomials 
    lapack_int n    = siNow.getSize();
    lapack_int nrhs = n;
    lapack_int lda  = n;
    lapack_int ldb  = nrhs;

    double * a = new double [n*n] ();
    double * b = new double [n*nrhs] ();
    lapack_int * p = new int [n] ();

    for (int cell=0; cell<n; cell++){
        // Cell indice  
        indice currentCell = start + siNow(cell);
        vector<vertex> work;

        for (auto & c: mi.faceCorner){
            int sj = currentCell[1] + c[1] + mi.vertexGhostLayerSize;
            int si = currentCell[0] + c[0] + mi.vertexGhostLayerSize;

            int fulllocali = si - mi.MPIlocalCellStart[0];
            int fulllocalj = sj - mi.MPIlocalCellStart[1];

            work.push_back(mi.lmesh[fulllocalj*mi.MPIlocalVertexSizeFull.at(0)+fulllocali]);
        }

        for (int r=0; r<n; r++){
            int xpow = r%siNow.getSize(0);
            int ypow = r/siNow.getSize(1);
            a[cell*n + r] = NumIntegralFace(work, {xpow, ypow}, center, scale, basePoly);
        }
    }

    fill(b,b+n*nrhs,0);
    for (int i=0; i<nrhs; i++) {b[i*n+i]=a[n*i];}

    int err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, p, b, ldb);

    if (err){
        printf("ERROR: Weno Basis Coefficient for order %d, %d. Error type %d \n",
                 siNow.getSize(0),siNow.getSize(1),err);
    }

    tensorpoly.setSize({siNow.getSize(0), siNow.getSize(1)});

    for (int p=0; p<n; p++){
        double * tmpcoef = new double [n]();
       
        for (int r=0; r<n; r++){
            tmpcoef[r] = b[r*n+p];
        }

        tensorpoly(p) = polynomial(siNow.getSize(0), siNow.getSize(1));
        tensorpoly(p).setCoef(tmpcoef);

        delete [] tmpcoef;
    }

    delete [] a;
    delete [] b;
    delete [] p;

    return 1;
}
