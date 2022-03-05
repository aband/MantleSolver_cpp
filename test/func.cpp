#include "func.h"

/*
 *Define functions used in calculation here
 */
double funcF(double t, valarray<double> target, double u){

    return u*u/2;
}

double funcG(double t, valarray<double> target, double u){

    return u*u/2;
}

double dfuncF(double t, valarray<double> target, double u){

    return u;
}

double dfuncG(double t, valarray<double> target, double u){

    return u;
}

double Initial_Condition(valarray<double>& target, const vector<double>& param){

    if (target[0]<0.5 && target[1]<0.5){
        return 0.5;
    } else if (target[0]>0.5 && target[1]<0.5){
        return 0.8;
    } else if (target[0]<0.5 && target[1]>0.5){
        return -0.2;
    } else {
        return -1.0;
    }
}
