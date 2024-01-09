#ifndef error_H_
#define error_H_

// Return error measured in energy norm or any arbitrary norm

PetscErrorCode L2ErrorElem(const vector<double>& coeff,
                           const bndryVal& bndryval,
                           const indice& globalElemIndic,
                           double (*func)(const vertex& point, 
                                          const vector<double>& param),
                           const valarray<double>& gwf,
                           const vector<vertex>& gpf,
                           basis& basis_,
                           Hdivmixed& hdiv_);

#endif
