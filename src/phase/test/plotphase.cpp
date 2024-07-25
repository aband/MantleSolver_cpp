#include <iostream>
#include <petsc.h>

#include "eutectic.h"

using namespace EUTECTIC;

int main(int argc, char **argv){

    phaseState* pPtr = new phaseState();

    int state = pPtr->EvalPhaseRegion(0.5,2);

    // Write phase diagram out to file 
    

    return 0;
}
