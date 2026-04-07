#include "../lagrange_tmp.h"
#include <iostream>
#include <valarray>
#include <vector>

static double testFunc(std::valarray<double>& point){

    //return point[0]*point[0];
	 return sin(point[1]);
}

static double testDeriv(std::valarray<double>& point){

    return cos(point[1]);
}

int main(int argc, char ** argv){

    int degree = 3 + 1;

    const int halfPts = std::ceil((degree+1)/2.0);
    const int numPts =  halfPts * 2;

    cout << "The number of sampling points is " << numPts << endl;

    // Random starting point
	 std::valarray<double> middle {0.04, -0.1};

    // Sampling interval
    double dx = 0.5; 

    // Unitnormal vector (direction)
	 //std::valarray<double> unitnormal {0.3*sqrt(2), 0.4*sqrt(2)}; 

    //std::valarray<double> unitNormal {1.0, 0.0}; 
    std::valarray<double> unitNormal {0.0, 1.0}; 

    // Initialize 
    LagrangeBasisDeriv lagDer(numPts-1);

	 std::vector<double> samples;
    samples.resize(numPts);

    std::vector<std::valarray<double>> sampleP;
    sampleP.resize(numPts);

    for (int s=0; s<halfPts; s++){
        
        std::valarray<double> point0 = middle - (halfPts - 0.5 - s)*dx*unitNormal; 
	 	  std::valarray<double> point1 = middle + (0.5+s)*dx*unitNormal;

        sampleP.at(s) = point0;
        sampleP.at(halfPts + s) = point1;

	     samples.at(s) = testFunc(point0);
        samples.at(halfPts + s) = testFunc(point1);
 }

    for (const auto& p: sampleP){
        cout << p[0] << "  " << p[1] << endl ;
	 }

    double work = 0.0;

    for (int i=0; i<numPts; i++){
        work -= lagDer.middle(numPts-1, i)/dx * samples.at(i);
	 }

    cout << "Analytical solution of the derivative is : " << testDeriv(middle) << endl;

    cout << "Numerical solution of the derivative is : " << work << endl;

    cout << "Standarf error funciton " << std::erf(0.5) << endl;

    return 1;
} 
