#include "diffusiveflux.h"
#include "trans_param.h"
#include "lagrange_tmp.h"

// Initial consition
double func(const vertex& point,
            const vector<double>& param){

    // Initial condition

    // Initialize with simple Reimann shock and rarefaction function
    // time inputed as param[0] 

    // Rarefraction initial condition
    //return 0.1*(3-point[0]);

    return 0.1;
}

// Burgers for testing
double diffunc(const double& u){

    // A Burgers type flux

    return u;
}

int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts){

    Tensor_zero(vertedge);
    Tensor_zero(horiedge);

    for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1]+1; j++){



    }

    return 1;
}
