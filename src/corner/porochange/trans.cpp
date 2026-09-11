#include "couple.h"
#include "transport.h"
// Transport boundary  condition for thermal and compositional conditions

// Initialize dimensionless composition and enthalpy
double InitCD(const valarray<double>& point,
              const vector<double>& param){

    // Constant composition value 
    //return 0.00;
	 //return 1.0; 
	 //if (point[0] >= 0.2 && point[0] <= 0.3){
    //    return 1.0;
	 //} else {
    //    return 0.0;
	 //}

    //return 1.0;
	 return 3.35;
    //return 3.0 - point[1]*0.12;
}

double InitHD(const valarray<double>& point,
              const vector<double>& param){

    // Linear simple distribution of enthalpy
	 // We pass nondimensionalize normalization factor in param.at(0)
    //double HD = 0.01;

    //double HD = 2.9-2.5*point[1];

    //if (point[1] < -0.20){HD = 2.9 + 2.5*0.20;}

    double HD = 0.01 - point[1]*0.12;

//    return 0.1;

    return HD;
}

double dfdu(double fneg, double fpos, double uneg, double upos, double alpha){
    // linear transport
    return 1.0;
}

double advfunc(const double& u, const vertex& vel, const vertex& unitnormal){

    return u*(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
}

double dfdu(const double& u){

    return 1.0;
}

// ====================================================================================
// Boundary condition for transport problem

// Check if it is outflow first then assign boundary condition to this transport problem
// This is advection boundary
int diriBndry(const vector<vertex>& points, 
                    vector<double>& value,
                    int flag){

    // flag switch between different variables.

    for (int g=0; g<points.size(); g++){

        if (flag == 0){

            //if (points.at(g)[1] == -0.5 && points.at(g)[0]>0.2 && points.at(g)[0]<0.3){
            if (points.at(g)[1] == -0.5){
                //value.at(g) = 1.0;
                value.at(g) = 3.35;

            } else {
                value.at(g) = 0.0;
            }

        } else {

            value.at(g) = 0.01;

        }

    }

    return 1;
}

// Diffusion boundary 
double diffBndry(const vertex& points, const vector<double>& param){

    // A cooling surface diffusion boundary condition

    if (points[1] > -0.0001){
  
        return 0.0;
        //return -2.3; // Zero celcius

	 } else {

        //return 10 + points[1]* param.at(0);
		  return 0.0;
	 }

}

int diffBndryType(const vertex& points){

    // 0: Zero diffusion boundary
    // 1: Dirichlet diffusion boundary

    if (points[1] > -0.0001){
        return 1;
		  //return 0; 
    } else {
        return 0;
    }

}
