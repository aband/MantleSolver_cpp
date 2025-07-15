#include "stencilpolynomial.h"

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
             tensorsigma({p,p}) = polyintegral(corners, area, center, h, p);
        }

    }

    return 1;
}
