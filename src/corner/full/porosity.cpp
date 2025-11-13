#include "couple.h"

static double AssignPorosity(const vertex& point){

    // used to identify incorrect porosity

    if (point[1] > -0.05){
        return 0.1;
    }else if (point[1] > -0.15 && point[1] < -0.05){
        return 0.05;
    }else if (point[1] > -0.25 && point[1] < -0.15){
        return 0.05;
    } else {
        return 0.0;
    }

}

int couple::computePorosity(){

    // Compute porosity on each gauss points 
    // With a fixed function
    edgeporo.clear();
	 int total = edgegauss.size();
    edgeporo.resize(total);

    for (int g=0; g<total; g++){
        edgeporo.at(g) = AssignPorosity(edgegauss.at(g)); 
    }

    return 1;
}

int couple::computePorosity_phase(){

    


    return 1;
}
