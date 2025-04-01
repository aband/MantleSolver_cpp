#include "driver.h"

// Read output HD and CD data form 
// output file as a restart
int ReadValues(const char * fieldname, int mark, vector<double>& data){

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    FILE * sol = fopen(filename, "r");

    for (int i=0; i<data.size(); i++){
        double val = 0.0;
        fscanf(sol, "%lf ", &val);
        data.at(i) = val;
        printf("%.16f \n", val);
    }

    fclose(sol);

    return 1;
}

int ReadValues(const char * filename, vector<double>& data){

    FILE * sol = fopen(filename, "r");

    for (int i=0; i<data.size(); i++){
        double val = 0.0;
        fscanf(sol, "%lf ", &val);
        data.at(i) = val;
        printf("%.16f \n", val);
    }

    fclose(sol);

    return 1;
}
