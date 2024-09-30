#include <iostream>
#include <fstream>
#include <petsc.h>

#include "eutectic.h"

using namespace EUTECTIC;

// Define a function relating Enthalpy with depth
double HDz(double z){

   double HD = z/10000;// rescale depth

   HD = 0.18 + 0.000001*z;

   return HD;
}

int main(int argc, char **argv){

    phase* pPtr = new phase();

    FILE *gridCD = fopen("gridCD.dat", "w");
    FILE *gridHD = fopen("gridHD.dat", "w");
    FILE *TD     = fopen("TD.dat", "w");
    FILE *Vf   = fopen("Vf.dat", "w");

    int seed = 50;

    double lithoP  = 5.2356e-06/10e-7;

    double deep = 60000.0;

    double HD_i = -0.2 + 10e-7*lithoP*deep, HD_f=1.3 + 10e-7*lithoP*deep, CD_i=0, CD_f=1.0; 

    for (int i=0; i<seed+1; i++){
        double CD = CD_i + i*(CD_f-CD_i) /(double)(seed) ;
        for (int j=0; j<seed+1; j++){
            double HD = HD_i + j*(HD_f- HD_i)/(double)(seed);

            fprintf(gridCD, "%f ", CD);
            fprintf(gridHD, "%f ", HD);

            pPtr->evalPhase(HD, CD, lithoP*60000);

            fprintf(TD, "%f ", pPtr->TD);
            fprintf(Vf, "%f ", pPtr->phi.mlt);

        }

        fprintf(gridCD, "\n");
        fprintf(gridHD, "\n");
        fprintf(TD, "\n");

    }
   
    FILE *TDz    = fopen("TDz.dat", "w");
    FILE *Vfz    = fopen("Vfz.dat", "w");

    FILE *TDz3    = fopen("TDz3.dat", "w");
    FILE *Vfz3    = fopen("Vfz3.dat", "w");

    FILE *TDz4    = fopen("TDz4.dat", "w");
    FILE *Vfz4    = fopen("Vfz4.dat", "w");

    // Output of depth related data
    double h = deep/100.0;

    for (int k=0; k<100; k++){
        double HD = HDz(k*h); 
		  double P = lithoP*k*h;
        pPtr->evalPhase(HD, 0.2, P);
        //std::cout << pPtr->phaseSplit(HD, 0.5, lithoP*h*k) << "  " << HD << std::endl;
        fprintf(TDz, "%f ", pPtr->TD);
        fprintf(Vfz, "%f ", pPtr->phi.mlt);

        std::cout <<P << "  " << pPtr->phi.mlt << std::endl;

        pPtr->evalPhase(HD, 0.4, P);
        fprintf(TDz3, "%f ", pPtr->TD);
        fprintf(Vfz3, "%f ", pPtr->phi.mlt);

        pPtr->evalPhase(HD, 0.6, P);
        fprintf(TDz4, "%f ", pPtr->TD);
        fprintf(Vfz4, "%f ", pPtr->phi.mlt);

    }

    FILE *HDzgrid = fopen("HDzgrid.dat", "w");
    FILE *Pzgrid = fopen("Pzgrid.dat", "w");
    FILE *TDz2    = fopen("TDz2.dat", "w");
    FILE *Vfz2    = fopen("Vfz2.dat", "w");

    // H-P slice on C=0.5
    for (int m=0; m<100; m++){
        double HD = HDz(m*h);
        for (int n=0; n<100; n++){
            double P  = lithoP*n*h;
            pPtr->evalPhase(HD, 0.5, P);

            fprintf(HDzgrid, "%f ", HD);
            fprintf(Pzgrid, "%f ", P);

            fprintf(TDz2,"%f ", pPtr->TD);
            fprintf(Vfz2, "%f ", pPtr->phi.mlt);

        }

        fprintf(HDzgrid, "\n");
        fprintf(Pzgrid, "\n");
        fprintf(TDz2,"\n");
        fprintf(Vfz2, "\n");
 
    }

    fclose(gridCD);
    fclose(gridHD);
    fclose(TD);
    fclose(Vf);
    fclose(TDz);
    fclose(Vfz);
    fclose(HDzgrid);
    fclose(Pzgrid);
    fclose(TDz2);
    fclose(Vfz2);

    fclose(TDz3);
    fclose(Vfz3);
    fclose(TDz4);
    fclose(Vfz4);

    return 0;
}
