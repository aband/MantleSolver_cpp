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

int getQuadVelVert(vertexSet& quadvel, 
                   const vector<vertex>& velocityField, int i, int j, int M, int N){

    return 1;
}

int getQuadVelHori(vertexSet& quadvel, 
                   const vector<vertex>& velocityField, int i, int j, int M, int N){

    return 1;
}
