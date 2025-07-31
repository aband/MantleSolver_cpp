#include "tensorstencilpoly.h"
static int getcenter(const vector<vector<vertex>>& cornerSet,
                     vector<vertex>& refcell,
                     const int& sizex, const int& sizey,
                     vertex& center, const double& h){

    // Get index for element 3
    int index = (sizey-1)*sizex;

    center = (cornerSet.at(0)[0]             + 
              cornerSet.at(sizex-1)[1]       +
              cornerSet.at(sizex*sizey-1)[2] + 
              cornerSet.at(index)[3] )/4.0   ;

    // Create reference cell corners
    refcell.clear();
    refcell.resize(4);
    vertex add = {-h/2, -h/2};
    refcell[0] = center + add;
    add = {h/2, -h/2};
    refcell[1] = center + add;
    add = {h/2, h/2};
    refcell[2] = center + add;
    add = {-h/2, h/2};
    refcell[3] = center + add; 

    return 1;
}

// Horner's method
static double horner(double x, const double* coef, int degree) {
  if(abs(x) <= 1) {

    double val = coef[degree];
    for(int i = degree-1; i >= 0; i--) {
      val = val*x + coef[i];
    }
    return val;

  } else {
   
    double val = coef[0];
    for(int i = 1; i <= degree; i++) {
      val = val/x + coef[i];
    }
    return pow(x,degree) * val;
    
  }
}

// ==========================================================================

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

tensorstencilpoly::~tensorstencilpoly(){

    if (coef) delete [] coef;

}

int tensorstencilpoly::setCoef(const MeshInfo& mi, const int& gstartx, const int& gstarty){

    // Extract tensor product stencil
    vector<vector<vertex>> cornerSet;

    for (int j=0; j<sizey; j++){
        for (int i=0; i<sizex; i++){
            indice global {i+gstartx, j+gstarty};
            cornerSet.push_back(extractCorners(mi, global));
        }
    }

    // Define 
    refarea = mi.L*mi.H/(double)(mi.MPIglobalCellSize[0]*
                                 mi.MPIglobalCellSize[1]);

    h = sqrt(refarea);

    getcenter(cornerSet, refcell, sizex, sizey, center, h); 

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
                //a[cell*n + j*sizex+i] = NumIntegralFace(work, {i,j}, center, h, basePoly);
                a[cell + (j*sizex+i)*n] = NumIntegralFace(work, {i,j}, center, h, basePoly);
            }
        }
    }

    fill(coef,coef+n*nrhs,0);
    for (int i=0; i<nrhs; i++) {coef[i*n+i]=a[i];}

    int err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, p, coef, ldb);

    if (err){
        printf("ERROR: Weno Basis Coefficient for order %d, %d. Error type %d \n",
                 sizex,sizey,err);
    }

    delete [] a;
    delete [] p;

    return 1;
}

double tensorstencilpoly::eval(const double& x,  const double& y, const int& ncell) const{

    double xx = (x-center[0])/h;
    double yy = (y-center[1])/h;

    double ycoef[sizey];

    int start = ncell*sizex*sizey;

    for (int j=0; j<sizey; j++){
        ycoef[j] = horner(xx, &coef[start], sizex-1);
        start += sizex;
    }

    return horner(yy, ycoef, sizey-1);
}

double tensorstencilpoly::eval(const double& x,  const double& y, 
                               const double& x0, const double& y0, 
                               const double& h,  const int& ncell) const{

    double xx = (x-x0)/h;
    double yy = (y-y0)/h;

    double ycoef[sizey];

    int start = ncell*sizex*sizey;

    for (int j=0; j<sizey; j++){
        ycoef[j] = horner(xx, &coef[start], sizex-1);
        start += sizex;
    }

    return horner(yy, ycoef, sizey-1);
}

int tensorstencilpoly::printCoef(){

    int n = sizex*sizey;

    for (int p=0; p<n; p++){
//        for (int r=0; r<n; r++){
//            cout << coef[r*n+p] << "  " ;
//        } cout << endl;
        printCoef(&coef[p*n], n);
		  cout << endl;
    } 

    return 1;
}

int tensorstencilpoly::printCoef(double * c, int n){

    for (int r=0; r<n; r++){
        cout << c[r] << "   ";
    }

    return 1;
}
