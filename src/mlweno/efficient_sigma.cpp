#include "stencilpolynomial.h"

double stencilpolynomial::polyintegral(const vector<vertex>& corners, 
                                       const double& area,
                                       const vertex& center, const double& h,
                                       const int& index1){

    double work = 0.0;

    const valarray<double>& gwf = GaussWeightsFace;
	 const vector<vertex>&   gpf = GaussPointsFace;

    vector<vertex> tmp = corners;           // Extract four corners
    for (auto & p: tmp) {p-=center; p/=h;}  // Transform points locally

    for (int i=0; i<gpf.size(); i++){
    
    }

    return work;
}

int stencilpolynomial::newsigmaintegral(const <vertex>& corners,
                                        const double& area,
                                        const vertex& center,
                                        const double& h,
                                        const int& max,
                                        Tensor<double>& deriv){

    

    return 1; 
}

// new efficient implementation of sigma
int stencilpolynomial::sigma(const vector<vertex>& corners, const double& area,
                             const vertex& center, const double& h, const string& type){

    // Assign sigma type
    sigmaType = type;

    // Compute basis tensor for Jiang-Shu sigma and polynomial sigma
    int total = 0;

    if (type == "JS"){
        // Jiang-Shu basis tensor
        total = size[0]*size[1];
        tensorsigma.setSize({total, total});
        for (int d=0; d<total; d++){
            tensorsigma({d,d}) = sigmaintegral(corners, area, center, h, i, j);}

        for (int j=0; j<total; j++){
            for (int i=0; i<j; i++){
                tensorsigma({i,j}) = sigmaintegral(corners, area, center, h, i, j);
                tensorsigma({j,i}) = tensorsigma({i,j});
            }
        }

    } else if (type == "Poly"){
        // Polynomial approximation of smoothness indicator sigma
        // Only the diagonal part is necessary
        total = ;
        tensorsigma.setSize(total, 1);

        for (int p=0; p<total; p++){
             tensorsigma({p,1}) = polyintegral(corners, area, center, h, p);
        }

    }

    return 1;
}

int 
