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

std::array<double, 3> trueSol(const vertex& point){

    array<double, 3> work;

    // ======================================================
    // Darcy test problem
    // return a predefined true solution
    // return <ux, uy, p> in this order
    // A manufactured solution satisfying Darcy equation

    // First scenerio
    // u + grad(p) = 0
    // div(u)      = 0
    // ux = -x/(x^2+y^2)
    // uy = -y/(x^2+y^2)
    // p  = 1/2 ln(x^2+y^2)

    //work[0] = -point[0]/(point[0]*point[0] + point[1]*point[1]);
    //work[1] = -point[1]/(point[0]*point[0] + point[1]*point[1]);
    //work[2] = 0.5*log(point[0]*point[0] + point[1]*point[1]);

    // Second scenerio
    // Divergence free linear velocity with arbitrary defined pressure field

    //work[0] = 1;   
    //work[1] = 1;
    //work[2] = point[0] * point[1];

    // Third scenerio
    work[0] = pow(point[0],2)*point[1];
    work[1] = -pow(point[1],2)*point[0];
    //work[2] = -point[0] + point[1];
    //work[0] = -point[0]*point[1];
    //work[1] = 0.5*pow(point[1],2);

    //work[2] = -0.5*point[0]*point[0] + 0.5*point[1]*point[1];
    work[2] = 1.0;

    // =================================================================
    // Test for Stokes problem
    //work[0] = cos(point[0])*sin(point[1]);
    //work[1] = -sin(point[0])*cos(point[1]);

    //work[2] = sin(point[0])*sin(point[1]);

    // Constant true solution
    //work[0] = pow(point[1],3);
    //work[1] = 0.0;
    //work[0] = 1;
    //work[1] = 1;
    //work[2] = 0.0;

    return work;
}

const vertex darcyPressureGrad(const vertex& point){

    // Auxiliary function.
    // Returns the gradient of scalar pressure field
    return {0.0,0.0};
}

const vertex stokesPressureGrad(const vertex& point){

    //return {cos(point[0])*sin(point[1]),
    //        sin(point[0])*cos(point[1])};

    return {0.0,0.0};
}

const vertex divdivVel(const vertex& point){

    //return {-2*cos(point[0])*sin(point[1]),
    //         2*sin(point[0])*sin(point[1])};

    //return {2*point[1],-2*point[0]};
    return {6*point[1],0.0};
}

const vertex Dirichlet_val(const vertex& point){

    std::array<double, 3> truesol = trueSol(point);

    // Test of Darcy part
    return {truesol[0], truesol[1]};
}

const vertex darcyForce(const vertex& point){

    // Return the arbitrarily defined right hand side
    // source term.
    std::array<double, 3> truesol = trueSol(point);

    vertex gradpressure = darcyPressureGrad(point);

    return {truesol[0] + gradpressure[0], truesol[1] + gradpressure[1]};
}

const vertex stokesForce(const vertex& point){

    return -1*divdivVel(point)+stokesPressureGrad(point);
}
