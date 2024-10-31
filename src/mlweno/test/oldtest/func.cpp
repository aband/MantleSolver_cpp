#include "func.h"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286

/*
 *Define functions used in calculation here
 */
double funcX(valarray<double>& target, const vector<double>& param){
    return 0.5*param[0]*param[0];
}

double funcY(valarray<double>& target, const vector<double>& param){
    return 0.5*param[0]*param[0];
}

double dfuncX(valarray<double>& target, const vector<double>& param){
    return param[0];
}

double dfuncY(valarray<double>& target, const vector<double>& param){
    return param[0];
}

double Initial_Condition(valarray<double>& target, const vector<double>& param){

	 // Oblique problem

	 if (target[0]<0.5 && target[1]<0.5){
		  return 0.5;
	 } else if (target[0]>0.5 && target[1]<0.5){
		  return 0.8;
	 } else if (target[0]<0.5 && target[1]>0.5){
		  return -0.2;
	 } else {
		  return -1.0;
	 }


/*
 *    // Diagnol testing
 *    if (target[0]+target[1]<1){
 *        return 0.5;
 *    }else {
 *        return 1.0;
 *    }
 *
 */


/*
 *    if (target[0]<0.0 && target[1]<0.0){
 *        return 0.5;
 *    } else if (target[0]<1.0 && target[0]>0.0 && target[1]<1.0 && target[1]>0.0){
 *        return 0.5;
 *    } else if (target[0]>1.0 && target[1]>1.0){
 *        return 0.5;
 *    } else if (target[0]>0.0 && target[1]<1.0){
 *        return -0.5;
 *    } else if (target[0]>1.0 && target[1]<1.0){
 *        return -0.5;
 *    } else {
 *        return 1.5;
 *    }
 *
 */

    //return pow(sin(PI*target[0]),2.0)+pow(sin(PI*target[1]),2.0);
    //return target[0];

// Riemann initial data
//    if (target[0]<0.5){
//        return 1.0;
//    }else {
//        return -0.5;
//    }
}



