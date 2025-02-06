#include "advectiveflux.h"
#include "trans_param.h"
#include "temp.h"

// Rarefraction case
// Initial consition
double func(const vertex& point,
            const vector<double>& param){

    // Initial condition

    // Initialize with simple Reimann shock and rarefaction function
    // time inputed as param[0] 

    // Rarefraction initial condition
    if (point[0] < 0.5 || point[0] >=(0.5*param[0]+1.5)){
        return 0;
    } else if (point[0]>=0.5 && point[0]<param[0]+0.5){
        return (point[0]-0.5)/param[0];
    } else if (point[0]>=point[0]+0.5 || point[0] <0.5*param[0]+1.5){
        return 1;
    } else {
        return 0;
    }
}

// Burgers for testing
double advfunc(const double& u, 
               const vertex& vel, const vertex& unitnormal){

    // A Burgers type flux

    return u*u/2.0 *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
}

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work){

    // compute df/du = df/dR * dR/du

    double direction = u*(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);

    unordered_map_arithmetic(work, du, std::plus<double>(), direction, std::multiplies<double>());

    return 1;
}

/**!
 * This function may be rewritten many times
 */
int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts){

    // Udpate every left and bottom edge for each cell
    //for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] + 1; j++){
    //for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] + 1; i++){
	 // Test for serial now
    Tensor_zero(vertedge);
    Tensor_zero(horiedge);
   
    const valarray<double>& gpe = GaussPointsEdge;
    vector<vertex> vel;
    vel.resize(gpe.size());
    for (int i=0; i<vel.size(); i++){
        vel.at(i) = {1,0};
    }

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double flux = 0.0;
        indice globalcell {i,j};
        indice cellout;

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, globalcell); 

        // Compute and restore Horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};
        // boundary
        if (j==0){
            // Temperatory
            flux    = 0.0;
        } else {
            cellout = globalcell + mi.faceNormal[0];
            flux    = edgefluxintegral(mi, globalcell, cellout, hori, allwgts, vel,
                                       ml, use, lu); 
        }

        horiedge({i,j}) = flux;
       
        vertexSet vert {corners.at(3), corners.at(0)};

        // boundary
        if (i==0){
            flux = 0.0;

        } else {
            cellout = globalcell + mi.faceNormal[3];
            flux    = edgefluxintegral(mi, globalcell, cellout, vert, allwgts, vel,
                                       ml, use, lu);
        }
        vertedge({i,j}) = flux;
    }}

    // Free outflow
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){

        indice gcell {mi.MPIglobalCellSize[0]-1, j};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet vert    = {corners.at(2), corners.at(1)};

        vertedge({mi.MPIglobalCellSize[0], j}) = edgefluxintegral(mi, gcell, vert, allwgts, vel, ml, use, lu);
    }

    return 1;
}

// Update edge flux and edge derivative of flux at once
int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   Tensor<derivative>& vertedgeder, Tensor<derivative>& horiedgeder,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts){

    const valarray<double>& gwe = GaussPointsEdge;
    vector<vertex> vel;
    vel.resize(gwe.size());
    for (int i=0; i<vel.size(); i++){
        vel.at(i) = {1,0};
    }

    // In this simplified test, no flux boundary are applied to all 
	 // four sides, flux = 0 and dflux = 0 as well
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double flux = 0.0;
        indice globalcell {i,j};
        indice cellout;

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, globalcell); 

        // Compute and restore Horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};
        // boundary
        if (j==0){
            // Temperatory
            flux    = 0.0;
        } else {
            cellout = globalcell + mi.faceNormal[0];
            edgefluxintegral(mi, globalcell, cellout, hori, 
            allwgts, vel, ml, use, lu, horiedgeder({i,j}), flux);
//            cout << i << "  " << j << endl;
//            unordered_map_print(horiedgeder({i,j})); cout << endl;
        }

        horiedge({i,j}) = flux;

        vertexSet vert {corners.at(3), corners.at(0)};

        // boundary
        if (i==0){
            flux    = 0.0;
        } else {
            cellout = globalcell + mi.faceNormal[3];
            edgefluxintegral(mi, globalcell, cellout, vert, allwgts, vel,
                             ml, use, lu, vertedgeder({i,j}), flux);

        }
        vertedge({i,j}) = flux;
    }}

    return 1;
}


