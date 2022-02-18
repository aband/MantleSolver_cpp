#include "weno.h"

void CreateSmoothnessIndicator(){

    for (auto & cell: index_set_stencil){
        sigmaS += pow(lsol[target_cell[1]][target_cell[0]] - lsol[cell[1]][cell[0]],2) ;
    }

    sigma *= 1.0/(double)(index_set_stencil.size()-1);

    sigma = pow(sigma, eta);
}
