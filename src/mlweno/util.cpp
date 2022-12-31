#include "util.h"

// Couple of constant functions
double constFunc(valarray<double>& point,const vector<double>& param){
    return 1.0;
}

double constFunc(valarray<double>& point,const vector<double>& param, double c){
    return c;
}

double constFunc(){
    return 1.0;
}

double constFunc(double c){
    return c;
}


// Evaluation of factorial
/*
 *int factorial(int top, int bottom){
 *    assert(top>bottom || top==bottom);
 *    if (top==bottom){ return 1;}
 *    else{ return top*factorial(top-1,bottom);}
 *}
 */

/*
 *int factorial(int top){
 *    assert(top>0 || top==0);
 *    if (top==0){ return 1;}
 *    else{ return top*factorial(top-1);}
 *}
 */

int factorial(int n){

    assert(n>0 || n==0);

    int * work = new int [n+1] ();

    work[0] = 1;
    for (int i=1; i<=n; i++){
        work[i] = i*work[i-1];
    }
    return work[n];

}

int factorial(int n, int m){
    assert(n>m || n==m);

    int * work = new int [n-m] ();

    work[0] = m+1;
    for (int i=1; i<n-m; i++){
        work[i] = (m+1+i)*work[i-1];
    }
    return work[n-m-1];
}

// Evalutaion of polynomial using Horner's method
double polyEval(double x, double * coef, int degree){

    if (abs(x) <= 1){

        double work = coef[degree];

        for (int r=degree-1; r>=0; r--){
            work = work*x + coef[r]; 
        } 

        return work;

    } else {

        double work = coef[0];

        for (int r=1; r<=degree; r++){
            work = work/x + coef[r];
        }

        return pow(x,degree)*work;
    }

}
