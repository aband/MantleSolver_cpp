#ifndef error_H_
#define error_H_

// Return error measured in energy norm or any arbitrary norm

PetscErrorCode L2ErrorElemInterior(const vector<double>& coeff,
                                   const indice& globalElemIndic,
                                   double (*func)(const vertex& point, 
                                                  const vector<double>& param),
                                   const valarray<double>& gwf,
                                   const vector<vertex>& gpf,
                                   basis& basis_,
                                   Hdivmixed& hdiv_);

PetscErrorCode L2ErrorElemBndry();

#endif
