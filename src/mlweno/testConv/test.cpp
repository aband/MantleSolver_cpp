#include <iostream>
#include <petsc.h>
#include "integral.h"

#include "stencil.h"
#include "util.h"
#include "polynomial.h"
//#include <adolc/adolc.h>

extern "C"{
#include "mesh.h"
#include "output.h"
}

using namespace std;

int main(int argc, char **argv){

    double coef[6] = {0.5,0.4,12,1,0.01,5};

    int maxDegree[2] = {2,3};

    //cout << polyEval(2.5, coef, 2) << endl;


    MLWENO::stencil <MLWENO::basisPolynomial *> stencilPoly(2,2);

    stencilPoly(1,1) = new MLWENO::basisPolynomial(maxDegree, coef);

    cout << stencilPoly(1,1)->eval(0.5,0.6) << endl;

    return 0;
}
