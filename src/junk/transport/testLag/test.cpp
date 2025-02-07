#include "lagrange_tmp.h"
#include <iostream>

int main(int argc, char **argv){

    LagrangeBasisDeriv lagDer(3);

    cout << lagDer.maxDegree() << endl;

    cout << lagDer.middle(3,1) << endl;

    cout << lagDer.middle(3,3) << endl;

    return 0;
}
