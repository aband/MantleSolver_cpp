#include "transport.h"

int TransportVariable::Print(const MeshInfo& mi, const char * filename){

    FILE * sol = fopen(filename,"w");

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        for (int e=0; e<4; e++){
        for (int g=0; g<3; g++){

            fprintf(sol, "%.12f ", my_recon.at(j*mi.MPIglobalCellSize[0] + i)->elem_val.at(e*3+g));

        }} fprintf(sol, "\n");

    }fprintf(sol,"\n");}

    fclose(sol);

    return 1;
}
