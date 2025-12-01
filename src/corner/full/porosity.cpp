#include "couple.h"

static double AssignPorosity(const vertex& point){

    // used to identify incorrect porosity

    return point[1];

}

int couple::computePorosity(){

    // Compute porosity on each gauss points 
    // With a fixed function
    edgeporo.clear();
    edgeporo.resize(edgegauss.size());

    for (int g=0; g<edgegauss.size(); g++){
        edgeporo.at(g) = AssignPorosity(edgegauss.at(g)); 
    }

    cellporo.clear();
    cellporo.resize(cellgauss.size());

    for (int g=0; g<cellgauss.size(); g++){
        cellporo.at(g) = AssignPorosity(cellgauss.at(g));
    }

    average_poro.clear();
    average_poro.resize(M_*N_);

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    double work = 0.0;
    double area = 0.0;

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        work = 0.0;
        area = 0.0;

        for (int g=0; g<gwf.size(); g++){

            double jac = 1;
        }

    }}

    return 1;
}

int couple::computePorosity_phase(){

    


    return 1;
}
