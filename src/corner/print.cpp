#include "couple.h"

int couple::printGaussPoints(){

    FILE * vertggridx = fopen("vertgaussgridx.dat", "w");
    FILE * vertggridy = fopen("vertgaussgridy.dat", "w");

    FILE * horiggridx = fopen("horigaussgridx.dat", "w");
    FILE * horiggridy = fopen("horigaussgridy.dat", "w");
//cout << edgegauss.size();
    // Vertical points first
    for (int j=0; j<N_; j++){
    for (int i=0; i<M_+1; i++){
      
        int dof = (j*(M_+1) + i)*3;
cout << dof << endl;
 
        for(int g=0; g<3; g++){
            //fprintf(vertggridx, "%e ", edgegauss.at(dof+g)[0]);
            //fprintf(vertggridy, "%e ", edgegauss.at(dof+g)[1]);
        }

    }fprintf(vertggridx, "\n ");
     fprintf(vertggridy, "\n ");}

    int tolvert = N_*(M_+1)*3;

    // Horizontal points second
    for (int j=0; j<N_+1; j++){
    for (int i=0; i<M_;   i++){

        int dof = tolvert + (j*M_ + i)*3;
cout << dof << endl;
        for (int g=0; g<3; g++){
            fprintf(horiggridx, "%e ", edgegauss.at(dof+g)[0]);
            fprintf(horiggridy, "%e ", edgegauss.at(dof+g)[1]);
        }

    }fprintf(horiggridx, "\n ");
     fprintf(horiggridy, "\n ");}

    fclose(vertggridx);
    fclose(vertggridy);
    fclose(horiggridx);
    fclose(horiggridy);

    return 1;
}
