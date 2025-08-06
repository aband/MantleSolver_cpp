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

// 1D Horner's method for polynomial derivative evaluation (all up to der) for p(x/h)
static int horner_der(int der, double* val, double x, double h, 
					        const double* coef, int degree) {
  double xx = x/h;

//  for (int t=0; t<degree+1; t++){
//  cout << coef[t] << "  " ;
//  }
//cout << endl;
  for(int i = 0; i<= der; i++) val[i] = 0;
  
  for(int i = degree; i >= 0; i--) {
    for(int d=der; d>=1; d--) val[d] = val[d]*xx + d*val[d-1];
    val[0] = val[0]*xx + coef[i];
  }

  for(int d=1; d<=der; d++) val[d] /= pow(h,d);

  return 1;
}

// 1D Horner's method for all polynomial derivative evaluation for p(x/h)
static int horner_der(double* val, double x, double h, 
					        const double* coef, int degree) {
  double xx = x/h;
  
  for(int i = 0; i<= degree; i++) val[i] = 0;
  
  for(int i = degree; i >= 0; i--) {
    for(int d=degree - i; d>=1; d--) val[d] = val[d]*xx + d*val[d-1];
    val[0] = val[0]*xx + coef[i];
  }

  for(int d=1; d<=degree; d++) val[d] /= pow(h,d);

  return 1;
}

// Evaluate val[m + (derX+1)*n] = D_x^m D_y^n p(x), p(x) = Sum_ij c_ij (x-x0)^i (y-y0)^j / h^(i+j)
static int polynomial2D_ders(int derX, int derY, double* val,
			                    double x, double y, double x0, double y0, double h,
			                    int polyn_degree, double* my_coef) {
  double xx0 = x-x0;
  double yy0 = y-y0;

  double valX[derX+1];
  double valY[derY+1];
  double yCoef[derX+1][polyn_degree+1];
  
  // Horner's method in x, for each power of y
  int sz = polyn_degree+1;
  int start = 0;
  for(int j = 0; j <= polyn_degree; j++) {
    horner_der(derX,valX,xx0,h,&my_coef[start],sz-1);
    for(int d = 0; d <= derX; d++) 
	 {yCoef[d][j] = valX[d];}
    start += sz;
  }

  // Horner's method in y
  for(int dX = 0; dX <= derX; dX++) {
    horner_der(derY,valY,yy0,h,yCoef[dX],polyn_degree);
    for(int dY = 0; dY <= derY; dY++) {
      val[dX + (derX+1)*dY] = valY[dY];
    }
  }
  return 1;
}

static int tensorpoly_ders(int derX, int derY, int sizex, int sizey, double* val, 
                           double x, double y, double x0, double y0, double h, double* my_coef){
//cout << x0 << "  " << y0 << "  " << h << endl;
//  for (int c=0; c<9; c++){
//      cout << my_coef[c] << "  ";
//  }

  double xx0 = x-x0;
  double yy0 = y-y0;

  double valX[derX+1] = {0};
  double valY[derY+1] = {0};
  double yCoef[derX+1][sizey];
  
  // Horner's method in x, for each power of y
  int sz = sizex;
  int start = 0;
  for(int j = 0; j < sizey; j++) {
    horner_der(derX,valX,xx0,h,&my_coef[start],sz-1);
    for(int d = 0; d <= derX; d++) 
	 {yCoef[d][j] = valX[d];}
    start += sz;
  }

  // Horner's method in y
  for(int dX = 0; dX <= derX; dX++) {
    horner_der(derY,valY,yy0,h,yCoef[dX],sizey-1);
    for(int dY = 0; dY <= derY; dY++) {
      val[dX + sizex*dY] = valY[dY];
    }
  }

    return 1;
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
    if (sigmabase) delete [] sigmabase;
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

int tensorstencilpoly::der(int derX, int derY, int ncell, double * dp, double x, double y){

    //polynomial2D_ders(derX, derY, dp, x, y, center[0], center[1], h, order, &coef[ncell*sizex*sizey]);
    tensorpoly_ders(derX, derY, sizex, sizey, dp, x, y, center[0], center[1], h, &coef[ncell*sizex*sizey]);

    return 1;
}

int tensorstencilpoly::der(int derX, int derY, int ncell, double * dp, double x, double y, double x0, double y0, double scale){

    //polynomial2D_ders(derX, derY, dp, x, y, x0, y0, scale, order, &coef[ncell*sizex*sizey]);
    tensorpoly_ders(derX, derY, sizex, sizey, dp, x, y, x0, y0, scale, &coef[ncell*sizex*sizey]);

    return 1;
}
/*
static double complete_sum(double * deri, double * derj, int sizex, int sizey){

    double work = 0.0;

    return work;
}
*/
static double tensor_sum(double * deri, double * derj, int total, 
                         double area, double h, double gw, double jac){

    double work = 0.0;

    for (int t=1; t<total; t++){
        work += deri[t] * derj[t] *gw*jac *pow(h*h,t); 
    }

    return work;
}

static int lowertri(double * sigmatensor, double * all, int total, 
				        double area, double h, double gw, double jac){

    for (int d=0; d<total; d++){
        sigmatensor[d*total+d] += tensor_sum(&all[d*total],
                                             &all[d*total],
                                  total,area, h, gw, jac);
    }

    for (int j=1; j<total; j++){
        for (int i=0; i<j; i++){
            sigmatensor[j*total+i] += tensor_sum(&all[j*total],
                                                 &all[i*total],
                                      total,area, h, gw, jac);   
            sigmatensor[i*total+j] = sigmatensor[j*total+i];
        }
    }

    return 1;
}

int tensorstencilpoly::setSigma(){

    // Compute sigma as a complete polynomial
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

//    vector<vertex> tmp = refcell;          // Extract four corners
//    for (auto & p: tmp) {p-=center; p/=h;} // Transform points locally

    int total = sizex*sizey;

    sigmabase = new double [total * total] ();
    memset(sigmabase, 0, total*total);

    double all[total * total] = {0};

    for (int g=0; g<gpf.size(); g++){    

        valarray<double> mapped = GaussMapPointsFace(gpf[g],refcell);
        double jac = abs(GaussJacobian(gpf[g],refcell));
        double gw = gwf[g];

        for (int ncell = 0; ncell<total; ncell++){
            der(sizex-1, sizey-1, ncell, &all[ncell*total], mapped[0], mapped[1]);
        }

        lowertri(sigmabase, all, total, refarea, h, gw, jac);
    }

    for (int t=0; t<total*total; t++){
        sigmabase[t] /= refarea;
    }

    return 1;
}

// Evaluation and sigma 
// simplified serial version
double tensorstencilpoly::eval(double ** localsol,
                               const int& startx, const int& starty,
                               const double& x, const double& y) const{

    double work = 0.0;

    int cell = 0;
    for (int j=0; j<sizey; j++){
    for (int i=0; i<sizex; i++){
        cell = j*sizex+i; 
        work += localsol[starty+j][startx+i] * eval(x,y,cell);
    }}

    return work;
}

double tensorstencilpoly::sigma(double ** localsol,
                                const int& startx, const int& starty){

    double work = 0.0;

    int total = sizex*sizey;
    int cell1=0;
    int cell2=0;
    for (int j1=0; j1<sizey; j1++){
    for (int i1=0; i1<sizex; i1++){
        cell1 = j1*sizex + i1;
        for (int j2=0; j2<sizey; j2++){
        for (int i2=0; i2<sizex; i2++){
            cell2 = j2*sizex + i2;
            work += localsol[starty+j1][startx+i1] * 
                    localsol[starty+j2][startx+i2] * 
                    sigmabase[cell1*total + cell2];
        }}
    }}

    return work;
}

// Print functions
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

int tensorstencilpoly::printSigmaBase(){

    int total = sizex*sizey;
    for (int j=0; j<total; j++){
    for (int i=0; i<total; i++){
        cout << sigmabase[j*total + i] << "  " ;
    }cout << endl;}

    return 1;
}
