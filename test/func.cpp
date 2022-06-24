#include "func.h"

/*
 *Define functions used in calculation here
 */
double funcX(valarray<double> target, const vector<double>& param){

    return param[0]*param[0]/2;
}

double funcY(valarray<double> target, const vector<double>& param){

    return param[0]*param[0]/2;
}

double dfuncX(valarray<double> target, const vector<double>& param){

    return param[0];
}

double dfuncY(valarray<double> target, const vector<double>& param){

    return param[0];
}

double Initial_Condition(valarray<double>& target, const vector<double>& param){

/*
 *    // Oblique problem
 *
 *    if (target[0]<0.5 && target[1]<0.5){
 *        return 0.5;
 *    } else if (target[0]>0.5 && target[1]<0.5){
 *        return 0.8;
 *    } else if (target[0]<0.5 && target[1]>0.5){
 *        return -0.2;
 *    } else {
 *        return -1.0;
 *    }
 *
 */

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
	  *if (target[0]<0.5 && target[1]>0.5){
	  *    return -1.0;
	  *} else if (target[0]+target[1]<1){
	  *    return -0.5;
	  *} else {
	  *    return 1.0;
	  *}
	  */

/*
 *    if (target[0]<0.5 && target[1]>0.5){
 *        return -1.0;
 *    } else if (target[0]+target[1]<1){
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

    return sin(10*(target[0]+target[1]));
}



