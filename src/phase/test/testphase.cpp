#include <iostream>
#include <petsc.h>

#include "phase.h"

int main(int argc, char **argv){

    // Test HCDdiagram
    EUTECTIC::evalPhase* pPtr = new EUTECTIC::evalPhase(0.1,0.2);  

    pPtr->L         = 4e5;
    pPtr->Xe        = 0.7;
    pPtr->multTemp1 = 1350;
    pPtr->eutecticTemp = 1227;

    // Set upeutectic physical properties



    return 0;
}
