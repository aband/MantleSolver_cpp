#include "error.h"

PetscErrorCode L2ErrorElem(const vector<double>& coeff, 
                           const bndryVal& bndryval, 
                           const indice& globalElemIndic,
                           array<double,3> (*func)(const vertex& point),
                           const valarray<double>& gwf,
                           const vector<vertex>& gpf,
                           basis& basis_,
                           Hdivmixed& hdiv_){

    // Calculate the L2 Error on the given element 

    for (int g=0; g<gwf.size(); g++){
        //Loop through gauess quadrature points
        vertex mapped = GaussMapPointsFace(gpg[g], basis_.corners());
        double jac = abs(gaussJacobian(gpf[g], basis_.corners()));
        double gw = gwf[g]
        elemError += gw*jac*;
    }

    return PETSC_SUCCESS;
}
