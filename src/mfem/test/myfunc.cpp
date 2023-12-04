#include "myFunc.h"

double AssignPorosity(const vertex& point, const double& l){

//    if (point[1] < 12000 && abs(point[0]) < point[1] + l){
//        return 0.05*pow(1.0-point[1]/120000,2) * (1-abs(point[0])/(l+point[1]));
//    } else {
//        return 0.0;
//    }

    return 1.0;

}

// ===================================================

const vertex Dirichlet_val(const vertex& point){

//    return {0.0,0.0};

    // Test of Darcy part
    return {-point[0]/(point[0]*point[0] + point[1]*point[1]),
            -point[1]/(point[0]*point[0] + point[1]*point[1])};

}
