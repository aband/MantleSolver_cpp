#include <iostream>
#include <fstream>
#include <petsc.h>

#include "eutectic.h"

using namespace EUTECTIC;

// Define a function relating Enthalpy with depth
double HDz(double z){

   double HD = z/10000;// rescale depth

   HD = 0.2 + 0.0000003*z;

   return HD;
}

int main(int argc, char **argv){

    phase* pPtr = new phase();

    FILE *gridCD = fopen("gridCD.dat", "w");
    FILE *gridHD = fopen("gridHD.dat", "w");
    FILE *TD     = fopen("TD.dat", "w");
    FILE *Vf   = fopen("Vf.dat", "w");

    int seed = 50;

    double HD_i = -0.2, HD_f=1.3, CD_i=0, CD_f=1.0; 

    for (int i=0; i<seed+1; i++){
        double CD = CD_i + i*(CD_f-CD_i) /(double)(seed) ;
        for (int j=0; j<seed+1; j++){
            double HD = HD_i + j*(HD_f- HD_i)/(double)(seed);

            fprintf(gridCD, "%f ", CD);
            fprintf(gridHD, "%f ", HD);

            pPtr->evalPhase(HD, CD, 0);

            fprintf(TD, "%f ", pPtr->TD);
            fprintf(Vf, "%f ", pPtr->phi.mlt);

        }

        fprintf(gridCD, "\n");
        fprintf(gridHD, "\n");
        fprintf(TD, "\n");

    }
   
    FILE *TDz    = fopen("TDz.dat", "w");
 
    // Output of depth related data
    double lithoP  = 5.2356e-06/10e-7;

    double h = 60000.0/100.0;

    for (int k=0; k<100; k++){
        double HD = HDz(k*h); 
        pPtr->evalPhase(HD, 0.5, lithoP*k*h);
        std::cout << pPtr->phaseSplit(HD, 0.5, lithoP*h*k) << "  " << HD << std::endl;
        fprintf(TDz, "%f ", pPtr->TD);
    }

    fclose(gridCD);
    fclose(gridHD);
    fclose(TD);
    fclose(Vf);
    fclose(TDz);
 
    return 0;
}
