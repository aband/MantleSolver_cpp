#include "couple.h"
#include "transport.h"
// Transport boundary  condition for thermal and compositional conditions

// Initialize dimensionless composition and enthalpy
double InitCD(const valarray<double>& point,
              const vector<double>& param){

    // Constant composition value 
    return 0.04;
}

double InitHD(const valarray<double>& point,
              const vector<double>& param){

    // Linear simple distribution of enthalpy
	 // We pass nondimensionalize normalization factor in param.at(0)
    //double HD = 0.01;

    double HD = 2.9-2.5*point[1];

    if (point[1] < -0.20){HD = 2.9 + 2.5*0.20;}

    return HD;
}

double dfdu(double fneg, double fpos, double uneg, double upos, double alpha){
    // linear transport
    return 1.0;
}
