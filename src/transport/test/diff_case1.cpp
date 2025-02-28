#include "diffusiveflux.h"
#include "trans_param.h"
#include "lagrange_tmp.h"

// Initial consition
double func(const vertex& point,
            const vector<double>& param){

    if (abs(point[0]-0.5)<0.1 && 
        abs(point[1]-0.5)<0.1){
        return 1.0;
    } else {
        return 0.0;
    }

}

// Burgers for testing
double diffunc(const double& u){

    // A Burgers type flux

    return u;
}

// Burgers for testing
double advfunc(const double& u, 
               const vertex& vel, const vertex& unitnormal){

    return u *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);

//      return u*u/2.0 *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
}

double dfdu(const double& u){

      return 1;
}

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work){

    double direction = dfdu(u)*(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);

    unordered_map_arithmetic(work, du, std::plus<double>(), direction, std::multiplies<double>());

    return 1;
}

int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts){

    Tensor_zero(vertedge);
    Tensor_zero(horiedge);

    // No parallel implemented here
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double flux = 0.0;
        indice gcell {i, j};
        indice gcellout;

        vertexSet corners = extractCorners(mi, gcell);

        // Compute flux on horizontal edges 
        vertexSet hori {corners.at(0), corners.at(1)};

        // boundary
        if (j==0){
            // No flow boundary for now
            flux = 0.0;
        } else {
            gcellout = gcell + mi.faceNormal[0];
            flux     = edgefluxintegral(mi, gcell, gcellout, hori, allwgts, ml, use, lu, "all");
        }

        horiedge({i,j}) = 0.1*flux;

flux = 0.0;

        // Compute flux on vertical edges
        vertexSet vert {corners.at(3), corners.at(0)};

        // boundary
        if (i==0){
            // No flow boundary for now
            flux = 0.0;
        } else {
            gcellout = gcell + mi.faceNormal[3];
            flux     = edgefluxintegral(mi, gcell, gcellout, vert, allwgts, ml, use, lu, "all");
        }

        vertedge({i,j}) = 0.1*flux;

    }}

    return 1;
}

int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   Tensor<derivative>& vertedgeder, Tensor<derivative>& horiedgeder,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts){


    return 1;
}
