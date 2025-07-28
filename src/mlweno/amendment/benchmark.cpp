#include <iostream> 
#include <cmath>
#include <cstdlib>

#include "../polynomial.h"

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
static void horner_der(int der, double* val, double x, double h, const double* coef, int degree) {
  double xx = x/h;
  
  for(int i = 0; i<= der; i++) val[i] = 0;
  
  for(int i = degree; i >= 0; i--) {
    for(int d=der; d>=1; d--) val[d] = val[d]*xx + d*val[d-1];
    val[0] = val[0]*xx + coef[i];
  }

  for(int d=1; d<=der; d++) val[d] /= pow(h,d);
}

// 1D Horner's method for all polynomial derivative evaluation for p(x/h)
static void horner_der(double* val, double x, double h, const double* coef, int degree) {
  double xx = x/h;
  
  for(int i = 0; i<= degree; i++) val[i] = 0;
  
  for(int i = degree; i >= 0; i--) {
    for(int d=degree - i; d>=1; d--) val[d] = val[d]*xx + d*val[d-1];
    val[0] = val[0]*xx + coef[i];
  }

  for(int d=1; d<=degree; d++) val[d] /= pow(h,d);
}

// Evaluate val[m + (derX+1)*n] = D_x^m D_y^n p(x), p(x) = Sum_ij c_ij (x-x0)^i (y-y0)^j / h^(i+j)
void polynomial2D_ders(int derX, int derY, double* val,
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
    for(int d = 0; d <= derX; d++) yCoef[d][j] = valX[d];
    start += sz;
    sz--;
  }

  // Horner's method in y
  for(int dX = 0; dX <= derX; dX++) {
    horner_der(derY,valY,yy0,h,yCoef[dX],polyn_degree);
    for(int dY = 0; dY <= derY; dY++) {
      val[dX + (derX+1)*dY] = valY[dY];
    }
  }
}

int main(){

    // A n*n tensor product polynomial
    double my_coef[5*5];

    for (int i=0; i<25; i++){
        my_coef[i] = (double)i;
    }

    double val[5*5];

    polynomial2D_ders(3,3,val, 1.0,1.0,0.0,0.0,0.2,4,my_coef);

    // Benchmark code
    for (int i=0; i<25; i++){
        std::cout << val[i] << " " ;
		  if (i%5 == 4){
            std::cout << std::endl;
		  }
    }

    polynomial testp = polynomial(5,5);
    testp.setCoef(my_coef);
    testp.printCoef();

    return 1;
}
