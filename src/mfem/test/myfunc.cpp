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

array<double, 3> trueSol1(const vertex& point){

    // return a predefined true solution
    // return <ux, uy, p> in this order
    // A manufactured solution satisfying Darcy equation
    // u + grad(p) = 0
    // div(u)      = 0
    // ux = -x/(x^2+y^2)
    // uy = -y/(x^2+y^2)
    // p  = 1/2 ln(x^2+y^2)

    array<double, 3> work;

    work[0] = -point[0]/(point[0]*point[0] + point[1]*point[1]);
    work[1] = -point[1]/(point[0]*point[0] + point[1]*point[1]);
    work[2] = 0.5*log(point[0]*point[0] + point[1]*point[1]);

    return work;
}
