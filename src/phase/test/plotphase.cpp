#include <iostream>
#include <fstream>
#include <petsc.h>

#include "eutectic.h"

using namespace EUTECTIC;

int main(int argc, char **argv){

    phase* pPtr = new phase();

    FILE *gridCD = fopen("gridCD.dat", "w");
    FILE *gridHD = fopen("gridHD.dat", "w");

    FILE *TD     = fopen("TD.dat", "w");

    int seed = 50;

    double HD_i = -0.2, HD_f=1.3, CD_i=0, CD_f=1.0; 

    for (int i=0; i<seed+1; i++){
        double CD = CD_i + i*(CD_f-CD_i) /(double)(seed) ;
        for (int j=0; j<seed+1; j++){
            double HD = HD_i + j*(HD_f- HD_i)/(double)(seed);

            fprintf(gridCD, "%f ", CD);
            fprintf(gridHD, "%f ", HD);

            pPtr->evalPhase(HD, CD);

            fprintf(TD, "%f ", pPtr->TD);
        }

        fprintf(gridCD, "\n");
        fprintf(gridHD, "\n");
        fprintf(TD, "\n");

    }

    fclose(gridCD);
    fclose(gridHD);
    fclose(TD);

    return 0;
}
