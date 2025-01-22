#include "polynomial.h"
#include "tensor.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

int main(int argc, char ** argv){

    polynomial * testpoly = new polynomial(2,2);

    double testcoef[4] = {1,2,3,4};

    testpoly->setCoef(testcoef);
    testpoly->printCoef();

    testpoly->resetDegree(1,4);

    cout << testpoly->eval(1,0) << endl;
    cout << testpoly->eval(0,1) << endl;

    testpoly->evalDer(4,0,4,1,1,0);

    Tensor<double> mytensor = Tensor<double>(2);

    mytensor.setSize({3,2});

    return 0;
}
