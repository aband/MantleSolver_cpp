#include "polynomial.h"
#include "stencilpolynomial.h"

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

    int total = size[0] * size[1];
    tensorsigma.setSize({total, total});
    Tensor_zero(tensorsigma);

    vector<Tensor<double>> allder;
    allder.resize(total);

    Tensor<double> der = Tensor<double>(2);
    der.setSize(size);

    double sigmaijg = 0.0;

    for (int g=0; g<gpf.size(); g++){
        valarray<double> mapped = GaussMapPointsFace(gpf[g],tmp);
        double jac = abs(GaussJacobian(gpf[g],tmp));
        double gw = gwf[g];

        // Integration over a selected element
        for (int cell=0; cell<total; cell++){
            allder.at(cell).setSize(size);
            Tensor_zero(allder.at(cell));
            tensorpoly(cell).evalDer(size[0]-1, size[1]-1,
                                     mapped[0], mapped[1], 1.0, allder.at(cell));
        }

//        for (int jcell=0; jcell<total; jcell++){
//            for (int icell=0; icell<=jcell; icell++){
                // Compute derivative of a complete polynomial
                // instead of a tensor product one
//                sigmaij = 0.0; 
//                for (int j=0; j<order+1; j++){
//                    int iStart = (j==0) ? 1:0;
//                    for (int i=iStart; i<order+1-j; i++){
//                        sigmaij += pow(h,2*(i+j))/area * gw*jac* 
//                                   allder.at(jcell)({i,j}) * 
//                                   allder.at(icell)({i,j});
//                    }
//                }
//            }
//        }

        for (int j=0; j<total; j++){
            for (int i=0; i<j+1; i++){
                Tensor_multi_add(allder.at(j), allder.at(i), gw*jac, der);

                // sum der
                // For reference, using the wrong scaling factor
                sigmaijg = 0.0; 
                for (int e=1; e<total; e++){
                    sigmaijg += der(i) * pow(area/h*h,e);
                }

                // add to sigma base tensor
                // lower triangle
                tensorsigma({i,j}) += sigmaijg;
            }
        }
    }

    return 1;
}
