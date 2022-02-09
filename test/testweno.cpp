#include <iostream>
#include <petsc.h>
#include "../include/weno.h"
#include "../include/integral.h"
#include "../include/mesh.h"

#include <adolc/adolc.h>

using namespace std;

adouble power(adouble x, int n)
{
    adouble z = 1;
    if (n>0) {
        int nh = n/2;
        z = power(x,nh);
        z *= z;
        if (2*nh !=n) {z *= x;}
        return z;
    } else {
        if (n==0) {
            return z;
        } else {
            return 1.0/power(x,-n);
        }
    }
}

void testoutput(const valarray<double> test){
    valarray<double> get = test;
	 /*
     *cout << test.size() << endl;
     *cout << get.size() << endl;
	  */
	 for (auto &i: get){
		  i+=1;
	 }
	 /*
	  *for (auto &i: get){
	  *    cout << i << endl;
	  *}
	  */
}

inline valarray<double> testoutput2(valarray<double> test){
    valarray<double> gettest = test;
    for (auto &i: gettest){
        i = 7+i;
    }
    return gettest;
}

double func(valarray<double>& point){
    //return point[0]*point[0];
    return sin(point[0])+cos(point[1]);
}

int main(int argc, char **argv){

    PetscErrorCode ierr;
    PetscInitialize(&argc, &argv, NULL, NULL);

    DM dm;

    DMCreate(PETSC_COMM_WORLD,&dm);

    // valarray has dynamic memory allocation
    int k = 2;
    valarray<double> test(k);

    test[0] = 1.0;
    test[1] = 2.0;

    vector< valarray<int> > index { {-1,0,}, {-1,-1}, {0,1}, {1,0} };

    wenobasis * work;

    work = new wenobasis(test,index,2);

    testoutput(test);

    valarray<double> duplicate = testoutput2(test);
    duplicate = duplicate*duplicate;

    // Test for integral function
    vector< valarray<double> > corner {{0.0},{1,0.0},{1.5,1.5},{0.0,1}};
    valarray<double> center {0.0,0.0};
    double h = 1.0;
    
    double result = NumIntegralFace(corner,center,h,func);

    cout << result << endl;
	 /*
     *cout << duplicate[0] << endl;
     *cout << test[0] << endl;
     *cout << index[0][0] << endl;
	  */
/*
 *    cout << work[0].center[0] << work[0].center[1] << endl;
 *    for (auto & i: index){
 *        cout << i[0] << i[1] << endl;
 *    }
 *
 */
    return 0;
}
