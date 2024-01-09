#include "error.h"

PetscErrorCode L2ErrorElem(Vec * u, 
                           const bndryVal& bndryval, 
                           const indice& globalElemIndic){

    // Calculate the L2 Error on the given element 

    for (){
        //Loop through gauess quadrature points
        elemError += ;
    }

    return PETSC_SUCCESS;
}
