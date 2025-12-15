#include "couple.h"

static double AssignPorosity(const vertex& point){

    // used to identify incorrect porosity

    //return abs(point[1]);

    //return 0.2;

    double l0 = 316228;
    double l = 20/l0;

    if (abs(point[1]) < 120*1000/l0 && abs(point[0]) < abs(point[1]) + 20/l0){

        double value = 0.05*pow((120*1000/l0 - abs(point[1]))/(120*1000/l0),2) * 
                               (1-abs(point[0])/(abs(point[1])+l));

        return value;
    } else {
        return 0.0;
    }

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

        vector<vertex> corners = extractCorners(mi, {i,j});

        for (int g=0; g<gwf.size(); g++){

            vertex mapped = GaussMapPointsFace(gpf[g], corners);

            double jac = abs(GaussJacobian(gpf[g], corners));
            double gw = gwf[g];

            work += gw*jac*AssignPorosity(mapped);
            area += gw*jac;
        }

        average_poro.at(j*M_+i) = work/area; 

    }}

    return 1;
}

int couple::computePorosity_phase(){

    


    return 1;
}
