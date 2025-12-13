#include "couple.h"

int couple::printGaussPoints(){

    FILE * vertggridx = fopen("vertgaussgridx.dat", "w");
    FILE * vertggridy = fopen("vertgaussgridy.dat", "w");

    FILE * horiggridx = fopen("horigaussgridx.dat", "w");
    FILE * horiggridy = fopen("horigaussgridy.dat", "w");

    // Vertical points first
    for (int j=0; j<N_  ; j++){
    for (int i=0; i<M_+1; i++){
      
        int dof = (j*(M_+1) + i)*3;
 
        for(int g=0; g<3; g++){
            fprintf(vertggridx, "%e ", edgegauss.at(dof+g)[0]);
            fprintf(vertggridy, "%e ", edgegauss.at(dof+g)[1]);
				//cout <<dof+g << "  " << edgegauss.at(dof + g)[0] << "  " << edgegauss.at(dof + g )[1] << endl;
        }
    }}
//    }fprintf(vertggridx, "\n ");
//     fprintf(vertggridy, "\n ");}

    int tolvert = N_*(M_+1)*3;

    // Horizontal points second
    for (int j=0; j<N_+1; j++){
    for (int i=0; i<M_;   i++){

        int dof = tolvert + (j*M_ + i)*3;
        for (int g=0; g<3; g++){
            fprintf(horiggridx, "%e ", edgegauss.at(dof+g)[0]);
            fprintf(horiggridy, "%e ", edgegauss.at(dof+g)[1]);
        }
    }}

//    }fprintf(horiggridx, "\n ");
//     fprintf(horiggridy, "\n ");}

    fclose(vertggridx);
    fclose(vertggridy);
    fclose(horiggridx);
    fclose(horiggridy);

    return 1;
}

char * GetFilename(const char * fieldname, int mark){

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    return filename;
}

char * GetFilenameAdd(const char * fieldname, const char * add, int mark){

    char * filename = (char *)malloc(strlen(fieldname)+15+4);

    char n_char[15];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, add);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    return filename;
}

int couple::printedgeval(int mark, const vector<double>& val,
                                   const char * fieldname){

    FILE * file = fopen(GetFilename(fieldname, mark), "w");

    // Print scalar values on gauss quadrature points
	 // Vertical points first
    for (int j=0; j<N_  ; j++){
    for (int i=0; i<M_+1; i++){

        int dof = (j*(M_+1) + i)*3;

        for (int g=0; g<3; g++){
            fprintf(file,"%12f ", val.at(dof+g));
        }

    }}

    int tolvert = N_*(M_+1)*3;

    // Horizontal points next
    for (int j=0; j<N_+1; j++){
    for (int i=0; i<M_;   i++){

        int dof = tolvert + (j*M_ + i)*3;
        for (int g=0; g<3; g++){
            fprintf(file, "%12f ", val.at(dof+g));
        }
    }}

    fclose(file);
    return 1;
}

int couple::printcellval(int mark, const vector<double>& val,
                                   const char * fieldname){

    // Print values assigned to cell center
    FILE * file = fopen(GetFilename(fieldname, mark), "w");

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        fprintf(file, "%12f ", val.at(j*M_+i));

    }}

    fclose(file);
    return 1;
}

int couple::printedgeval(int mark, const vector<vertex>& val,
                                   const char * fieldname){

    FILE * filex = fopen(GetFilenameAdd(fieldname, "x", mark), "w");
    FILE * filey = fopen(GetFilenameAdd(fieldname, "y", mark), "w");

    // Print scalar values on gauss quadrature points
    // Vertical points first
    for (int j=0; j<N_  ; j++){
    for (int i=0; i<M_+1; i++){

        int dof = (j*(M_+1) + i)*3;

        for (int g=0; g<3; g++){
            fprintf(filex,"%12f ", val.at(dof+g)[0]);
            fprintf(filey,"%12f ", val.at(dof+g)[1]);
        }

    }}

    int tolvert = N_*(M_+1)*3;

    // Horizontal points next
    for (int j=0; j<N_+1; j++){
    for (int i=0; i<M_;   i++){

        int dof = tolvert + (j*M_ + i)*3;
        for (int g=0; g<3; g++){
            fprintf(filex, "%12f ", val.at(dof+g)[0]);
            fprintf(filey, "%12f ", val.at(dof+g)[1]);
        }
    }}

    fclose(filex);
    fclose(filey);

    return 1;
}

int couple::printedgeporosity(int mark){

    printedgeval(mark, edgeporo, "porosity");

    printcellval(mark, average_poro, "aveporo");

    return 1;
}

int couple::printedgevel(int mark, const vector<vertex>& stokesvel,
                                   const vector<vertex>& darcyvel){

    printedgeval(mark, stokesvel, "edgevel_stokes"); 

    printedgeval(mark, darcyvel , "edgevel_darcy");

    return 1;
}
