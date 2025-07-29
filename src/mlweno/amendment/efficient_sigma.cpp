#include "polynomial.h"

stencilpolynomial::stencilpolynomial(const int& order){
    size[0] = order+1;
    size[1] = order+1;

    tensorpoly.setSize(size);

    for (int i=0; i<tensorpoly.getSize(); i++){
        tensorpoly(i) = polynomial(size[0],size[1]);
    }
}

int stencilpolynomial::sigmacomplete(const vector<vertex>& corners,
					                      const double& area,
									          const vertex& center, 
									          const double& h){

    // Compute sigma as a complete polynomial
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    vector<vertex> tmp = corners;          // Extract four corners
    for (auto & p: tmp) {p-=center; p/=h;} // Transform points locally

    Tensor<double> der = Tensor<double>(2);
    der.setSize(size);
    Tensor_zero(der);

    int total = size[0] * size[1];
    tensorsigma.setSize({total, total});

    // Create two derivative tensor
    Tensor<double> der1 = Tensor<double>(2);
    der1.setSize(size);
    Tensor<double> der2 = Tensor<double>(2);
    der2.setSize(size);

    for (int g=0; g<gpf.size(); g++)
        valarray<double> mapped = GaussMapPointsFace(gpf[i],tmp);
        double jac = abs(GaussJacobian(gpf[i],tmp));
        double gw = gwf[i];

        for (int jcell=0; j<total; j++){
					 Tensor_zero(der1);
                tensorpoly(jcell).evalDer(size[0]-1, size[1]-1, 
                                          mapped[0], mapped[1], 1.0, der1);
            for (int icell=0; i<=j; i++){
                Tensor_zero(der2);
                tensorpoly(icell).evalDer(size[0]-1, size[1]-1, 
                                          mapped[0], mapped[1], 1.0, der2);

                // Compute derivative of a complete polynomial
					 // instead of a tensor product one
					 double sigmaij = 0.0; 
                for (int j=0; j<order+1; j++){
                    int iStart = (j==0) ? 1:0;
                    for (int i=iStart; i<order+1 -j; i++){
                        sigmaij += pow(area/(h*h),2*(i+j)) * gw*jac* 
										     der1({i,j}) * der2({i,j});
                    }	
                }
            }
        }
    }

    return 1;
}
