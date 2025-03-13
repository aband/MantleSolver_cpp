#include <iostream>
#include <fstream>
#include <petsc.h>

#include "eutectic_rescaled.h"

using namespace EUTECTIC;

// Define a function relating Enthalpy with depth
double HDz(double zD){

   double HD = 0.18;

   return HD;
}

int main(int argc, char **argv){

    phase* pPtr = new phase();

    FILE *gridCD = fopen("gridCD.dat", "w");
    FILE *gridHD = fopen("gridHD.dat", "w");

    int seed = 50;

    pPtr->printInfo();

    

    fclose(gridCD);
    fclose(gridHD);

    return 0;
}
