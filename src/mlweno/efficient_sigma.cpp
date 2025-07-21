#include "stencilpolynomial.h"

stencilpolynomial::stencilpolynomial(const int& inorder){

    // Create an even stencil with given order
    order = inorder;

    size.at(0) = order+1; 
    size.at(1) = order+1;

    tensorpoly.setSize(size);

    for (int i=0; i<tensorpoly.getSize(); i++){
        tensorpoly(i) = polynomial(size.at(0), size.at(1));
    }

	 // Define cover

}

stencilpolynomial::stencilpolynomial(const int& sizex,
                                     const int& sizey,
                                     const int& inorder){
    order = inorder;

    size.at(0) = sizex;
    size.at(1) = sizey;

    tensorpoly.setSize(size);

    for (int i=0; i<tensorpoly.getSize(); i++){
        tensorpoly(i) = polynomial(size.at(0), size.at(1));
    }
}

double stencilpolynomial::polyintegral(const vector<vertex>& corners, 
                                       const double& area,
                                       const vertex& center, const double& h,
                                       const int& index1){

    // Compute poly basis eta
    double work = 0.0;

    const valarray<double>& gwf = GaussWeightsFace;
	 const vector<vertex>&   gpf = GaussPointsFace;

    vector<vertex> tmp = corners;           // Extract four corners
    for (auto & p: tmp) {p-=center; p/=h;}  // Transform points locally

    for (int a=0; a<alpha; a++){
        for (int i=0; i<gpf.size(); i++){
            valarray<double> mapped = GaussMapPointsFace(gpf[i], tmp);
            double jac = abs(GaussJacobian(gpf[i], tmp));
            double gw  = gwf[i];

		      tensorpoly(a).evalDer(size.at(0)-1,size.at(1)-1,
							   	       mapped[0],mapped[1],1.0, der);
        }
    }

    return work;
}

int stencilpolynomial::newsigmaintegral(const <vertex>& corners,
                                        const double& area,
                                        const vertex& center,
                                        const double& h,
                                        const int& max,
                                        Tensor<double>& deriv){

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    int total = size.at(0)*size.at(1);
    derv.setSize({total, total});
    Tensor_zero(derv);

    vector<vertex> tmp = corners;          // Extract four corners
    for (auto & p: tmp) {p-=center; p/=h;} // Transform points locally

    int total = size[0] * size[1];  

    Tensor<double> der = Tensor<double>(2);
    der.setSize(size);
    Tensor_zero(der);

    // Not repeating calculation 
    for (int a=0; a<total; a++){	
        for (int i=0; i<gpf.size(); i++){

            valarray<double> mapped = GaussMapPointsFace(gpf[i], tmp);

            double jac = abs(GaussJacobian(gpf[i], tmp));
            double gw  = gwf[i];

		      tensorpoly(a).evalDer(size.at(0)-1,size.at(1)-1,
			    			            mapped[0],mapped[1],1.0, der);
        }
    }

    return 1; 
}

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

        for (int p=0; p<total; p++){tensorsigma({p,1}) = polyintegral(corners, area, center, h, p);}
    }
    return 1;
}
