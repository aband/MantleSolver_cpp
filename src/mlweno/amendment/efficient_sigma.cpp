#include "polynomial.h"

int stencilpolynomial::sigmaintegral(const vector<vertex>& corners,
					                      const double& area,
												 const vertex& center, const double& h,
												 vector<Tensor<double>>& allder){

    //Compute all at once
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    vector<vertex> tmp = corners;          // Extract four corners
    for (auto & p: tmp) {p-=center; p/=h;} // Transform points locally

    int total = size[0] * size[1];  

    Tensor<double> tmp = Tensor<double>(2);
    tmp.setSize(size);
    Tensor_zero(tmp);

    int alpha = size[0]*size[1];

    for (int i=0; i<gpf.size(); i++){
        valarray<double> mapped = GaussMapPointsFace(gpf[i],tmp);
        double jac = abs(GaussJacobian(gpf[i],tmp));
        double gw = gwf[i];

        for (int a=0; a<alpha; a++){
             tensorpoly(a).evalDer(size[0]-1,size[1]-1,
									        mapped[0], mapped[1], 1.0, tmp);
        }

    }

    return 1;
}
