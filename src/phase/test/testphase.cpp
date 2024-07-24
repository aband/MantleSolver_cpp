#include <iostream>
#include <petsc.h>

#include "phase.h"

int main(int argc, char **argv){

    // Test HCDdiagram
    EUTECTIC::evalPhase* pPtr = new EUTECTIC::evalPhase(6,0.2);  

    //pPtr->ViewPhysics();

    // Set upeutectic physical properties
    //pPtr->ViewPhase(); 

    //pPtr->EvalPhase(-0.5,0.2);

    //pPtr->ViewPhase();

    //pPtr->EvalPhase(-0.5,0.1);

    //pPtr->ViewPhase();

    //pPtr->EvalPhase(3,0.1);

    pPtr->ViewPhase();

    return 0;
}
