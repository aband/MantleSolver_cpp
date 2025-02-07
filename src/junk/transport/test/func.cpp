#include "func.h"

/**
 * Change functions for transport part here.
 * Define transport functions and derivatives
 * 2D Burger's equation
 *
 */
double funcX(vertex x, double u, double t){
    return 0.5*u*u;
}

double dfuncX(vertex x, double u, double t){
    return u;
}

double funcY(vertex x, double u, double t){
    return 0.5*u*u;
}

double dfuncY(vertex x, double u, double t){
    return u;
}

// ======== Diffusion =======================
double diffFunc(double u){
    return u;
}

double dDiffFunc(double u){
    return 1;
}

// ======== Initial distribution ============
double InitialDistribution(vertex& point, const vector<double>& param){

//    if (point[0] < 0){
//        return -1;
//    } else {
//        return 1;
//    }

	 //if (point[0]<-1.0/param[0]){
//		  return point[0]*point[0]+point[1]*point[1];
//	     return sin(point[0]*3.0)+cos(point[1]/2.0) + point[0]*(point[1]+1);
//	 } else {
//		  return point[0]*point[0]*point[1]*point[1] + 1.0;
//	     return sin(point[0]*3.0)+cos(point[1]/2.0) + point[0]*(point[1]+1) + 1;
//	 }

    //return sin(point[0]*3.0+0.5)+cos(point[1]/2.0-0.2) + pow(point[0]+0.1,3)*(point[1]+1);
    //return point[0]*point[0] + point[1]*point[1];
    //return point[0] + point[1];

    // Initial value for sine wave 2D Burger's equation
    return pow(sin(M_PI*(point[0]+1)/2),2)*pow(sin(M_PI*(point[1]+1)/2),2);

    //if (abs(point[0])+abs(point[1])<0.5){
    //    return 1;
    //} else {
    //    return 0;
    //}
}
