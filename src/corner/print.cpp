#include "couple.h"

int couple::printGaussPoints(){

    FILE * vertggridx = fopen("vertgaussgridx.dat", "w");
    FILE * vertggridy = fopen("vertgaussgridy.dat", "w");

    FILE * horiggridx = fopen("horigaussgridx.dat", "w");
    FILE * horiggridy = fopen("horigaussgridy.dat", "w");

    // Vertical points first
    for (int j=0; j<N_; j++){
    for (int i=0; i<M_+1; i++){
      
        int dof = (j*(M_+1) + i)*3;
        for(int g=0; g<3; g++){
            fprintf();
            fprintf();
        }

    }fprintf(vertggridx, "\n ");
     fprintf(vertggridy, "\n ");}

    // Horizontal points second
    for (int j=0; j<N_+1; j++){
    for (int i=0; i<M_;   i++){

        int dof = (j*M_ + i)*3;
        for (int g=0; g<3; g++){
            fprintf();
            fprintf();
        }

    }fprintf(horiggridx, "\n ");
     fprintf(horiggridy, "\n ")}

    fclose(vertggridx);
    fclose(vertggridy);
    fclose(horiggridx);
    fclose(horiggridy);

    return 1;
}
