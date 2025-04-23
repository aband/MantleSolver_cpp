#include <iostream>
#include <fstream>
#include <petsc.h>

#include "eutectic_rescaled.h"

using namespace EUTECTIC;

// Define a function relating Enthalpy with depth
double HDz(double zD){

   double HD = 3.0;

   return HD;
}

int main(int argc, char **argv){

    phase* pPtr = new phase();

    FILE *gridCD = fopen("gridCD.dat", "w");
    FILE *gridHD = fopen("gridHD.dat", "w");
    FILE *TD     = fopen("TD.dat", "w");
    FILE *Vf   = fopen("Vf.dat", "w");
    FILE *dTdH = fopen("dTdH.dat", "w");
    FILE *dTdC = fopen("dTdC.dat", "w");

    int seed = 100;

    pPtr->printInfo();

    double depth = 5;                //Dimensionless depth
    double h = depth / (double)seed; //Length interval

    double CD = 0.0;
    double HD = 0.0;

    double HD_i = pPtr->TDe0 - 0.2 , HD_f= pPtr->TDm0 + pPtr->LD;
    double CD_i=0, CD_f=pPtr->Xe; 
    // CD - HD plot
    for (int i=0; i<seed+1; i++){
        CD = CD_i + i*(CD_f-CD_i) /(double)(seed) ;
        for (int j=0; j<seed+1; j++){
            HD = HD_i + j*(HD_f- HD_i)/(double)(seed);

            // Plotting with shifted values
            fprintf(gridCD, "%f ", CD/pPtr->Xe);
            fprintf(gridHD, "%f ", HD-pPtr->TDe0);

            pPtr->evalPhase(HD, CD, 100000);

            fprintf(TD, "%f ", pPtr->pc.TDp - pPtr->TDe0);
            fprintf(Vf, "%f ", pPtr->pc.phil);

            fprintf(dTdC, "%f ", pPtr->pc.dTD_dCD);
            fprintf(dTdH, "%f ", pPtr->pc.dTD_dHD);

        }
        fprintf(gridCD, "\n");
        fprintf(gridHD, "\n");
        fprintf(TD, "\n");
        fprintf(dTdC, "\n");
        fprintf(dTdH, "\n");
    }

    // Depth plot

    fclose(gridCD);
    fclose(gridHD);
    fclose(TD);
    fclose(Vf);
    fclose(dTdH);
    fclose(dTdC);

    return 0;
}
