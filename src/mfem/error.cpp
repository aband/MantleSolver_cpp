#include "error.h"

double L2ErrorElemInterior(const vector<double>& weight, 
                           const indice& globalElemIndic,
                           std::array<double,3> (*func)(const vertex& point),
                           const valarray<double>& gwf,
                           const vector<vertex>& gpf,
                           basis& basis_,
                           Hdivmixed& hdiv_){

    // Calculate the L2 Error on the given interior element 

    for (int g=0; g<gwf.size(); g++){
        //Loop through gauess quadrature points
        vertex mapped = GaussMapPointsFace(gpg[g], basis_.corners());
        double jac = abs(gaussJacobian(gpf[g], basis_.corners()));
        double gw = gwf[g];
        // Evaluate all eight basis functions for the target element
		  std::array<vertex, 8> hdivwork = hdiv_.ComputeHdivmixed(basis_, mapped);
        // Combine these values with weights (calculated solution)
        valarray<double> approxVal = {0.0,0.0};
        for (int i=0; i<8; i++){
            approxVal += weight[i]*hdivwork[i]; 
        }
        // Get exact values
		  std::array<double,3> trueSol = func(mapped);

        valarray<double> diff = approxVal - trueSol;

        elemError += gw*jac*(diff[0]*diff[0] + diff[1]*diff[1]);
    }

    return pow(elemError,0.5);
}
