#include "driver.h"

int Driver::updateEdgeFluxDiff(Tensor<double>& vertedgeHD, Tensor<double>& horiedgeHD,
                               Tensor<double>& vertedgeCD, Tensor<double>& horiedgeCD,
                               const Tensor<weights>& allwgtsHD, double ** lHD,
                               const Tensor<weights>& allwgtsCD, double ** lCD){

    Tensor_zero(vertedgeHD);
    Tensor_zero(horiedgeHD);
    Tensor_zero(vertedgeCD);
    Tensor_zero(horiedgeCD);

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double flux = 0.0;

        indice gcell {i.j};
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

        horiedge({i,j}) = kdiff*flux;

        // =====================================================================================
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

        vertedge({i,j}) = kdiff*flux;

    }}

    return 1;
}
