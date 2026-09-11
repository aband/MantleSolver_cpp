#include "couple.h"

#include <random>

static double AssignPorosity(const vertex& point){

    double value = 0.0;

    std::uniform_real_distribution<double> real_distr(0.0, 0.3);

    std::random_device rd;
    std::mt19937 gen(rd());

    if (abs(point[1]) < 0.3 && abs(point[0]) < abs(point[1]) + 0.0005){

        value = real_distr(gen);

    } else {

        value = 0.0;

	 }
/*
    if (point[1] >-0.4 && point[1] < -0.3 && point[0]>0.4 && point[0]<0.51){

        value = 0.01;
	 }
*/
    return value;

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
