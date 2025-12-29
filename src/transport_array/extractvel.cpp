#include "extractvel.h"

int constVelField(vector<vertex>& vfield, int M, int N){

    int totalv = (N*(M+1) + (N+1)*M)*3;

    vfield.clear();

    vfield.resize(totalv);

    // Assign velocity across vertical edges first
    for (int j=0; j<N; j++){
    for (int i=0; i<M+1; i++){

        int prev = (j*(M+1) + i)*3;

        for (int g=0; g<3; g++){
            vfield.at(prev + g) = {1.0, 0.0};
        }

    }}

    // Assign velocity acorss horizontal edges next
    for (int j=0; j<N+1; j++){
    for (int i=0; i<M; i++){

        int prev = (N*(M+1) + j*M+i)*3;

        for (int g=0; g<3; g++){
            vfield.at(prev + g) = {1.0, 0.0};
        }

    }}

    return 1;
}

static char * obtainFilename(const char * fieldname, const char * tail, int mark){

    char * filename = (char *)malloc(strlen(fieldname)+20+4);

    char n_char[10];
    std::sprintf(n_char, "%d", mark);

    strcpy(filename, fieldname);
    strcat(filename, tail);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    return filename;
}

int printVelField(const vector<vertex>& velocityField, int M, int N, int mark, const char * fieldname){

    FILE * vertvx = fopen(obtainFilename(fieldname, "vertvx", mark),"w");
    FILE * vertvy = fopen(obtainFilename(fieldname, "vertvy", mark),"w");

    FILE * horivx = fopen(obtainFilename(fieldname, "horivx", mark),"w");
    FILE * horivy = fopen(obtainFilename(fieldname, "horivy", mark),"w");

    // Vertical velocity
    for (int j=0; j<N;   j++){
    for (int i=0; i<M+1; i++){

        int prev = (j*(M+1)+i)*3;

        for (int g=0; g<3; g++){
            fprintf(vertvx, "%e ", velocityField.at(prev + g)[0]); 
            fprintf(vertvy, "%e ", velocityField.at(prev + g)[1]);
        }
    }}

    // Horizontal velocity
    for (int j=0; j<N+1; j++){
    for (int i=0; i<M;   i++){

        int prev = (N*(M+1) + j*M + i)*3;

        for (int g=0; g<3; g++){
            fprintf(horivx, "%e ", velocityField.at(prev + g)[0]); 
            fprintf(horivy, "%e ", velocityField.at(prev + g)[1]);
        }
    }}

    fclose(vertvx);
    fclose(vertvy);
    fclose(horivx);
    fclose(horivy);

    return 1;
}

int printGaussGrid(int M, int N, const MeshInfo& mi){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> vertgaussp;
    vertgaussp.resize(gpe.size());

    std::vector<vertex> horigaussp;
    horigaussp.resize(gpe.size());

    vertexSet edge;

    FILE * vertggridx = fopen("vertgaussgridx.dat", "w");
    FILE * vertggridy = fopen("vertgaussgridy.dat", "w");

    FILE * horiggridx = fopen("horigaussgridx.dat", "w");
    FILE * horiggridy = fopen("horigaussgridy.dat", "w");

    vertexSet vertedge;
    vertexSet horiedge;

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        indice gcell {i,j};
        indice gcellout;
        vertexSet corners = extractCorners(mi, gcell);

        vertedge = {corners.at(0), corners.at(3)};
        horiedge = {corners.at(0), corners.at(1)};

        for (int g=0; g<gpe.size(); g++){
            vertgaussp.at(g) = GaussMapPointsEdge({gpe[g]},vertedge);
            horigaussp.at(g) = GaussMapPointsEdge({gpe[g]},horiedge);
        }   
 
        for (int g=0; g<gpe.size(); g++){

            fprintf(vertggridx, "%e ", vertgaussp.at(g)[0]);
            fprintf(vertggridy, "%e ", vertgaussp.at(g)[1]);

            fprintf(horiggridx, "%e ", horigaussp.at(g)[0]);
            fprintf(horiggridy, "%e ", horigaussp.at(g)[1]);
        }
    } fprintf(vertggridx, "\n ");
      fprintf(vertggridy, "\n ");
      fprintf(horiggridx, "\n ");
      fprintf(horiggridy, "\n ");}

    for (int i=0; i<M; i++){

        indice gcell {i, N-1};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet hori    = {corners.at(3), corners.at(2)};

        for (int g=0; g<gpe.size(); g++){
            horigaussp.at(g) = GaussMapPointsEdge({gpe[g]},horiedge);
        } 

        fprintf(horiggridx, "\n ");
        fprintf(horiggridy, "\n ");
    }

    for (int j=0; j<N; j++){

        indice gcell {M-1, j};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet vert    = {corners.at(1), corners.at(2)};

        for (int g=0; g<gpe.size(); g++){
            vertgaussp.at(g) = GaussMapPointsEdge({gpe[g]},horiedge);
        } 

        fprintf(vertggridx, "\n ");
        fprintf(vertggridy, "\n ");
    }

    fclose(vertggridx);
    fclose(vertggridy);
    fclose(horiggridx);
    fclose(horiggridy);

    return 1;
}

int getQuadVelVert(vertexSet& quadvel, 
                   const vector<vertex>& velocityField, 
						 int i, int j, int M, int N){

    return 1;
}

int getQuadVelHori(vertexSet& quadvel, 
                   const vector<vertex>& velocityField, 
						 int i, int j, int M, int N){

    return 1;
}
