#include "transfunc.h"

double dfdu(const double& u){

    return u*u*u;
}

// Burgers for testing
double advfunc(const double& u, 
               const vertex& vel, const vertex& unitnormal){

    // A Burgers type flux
//    return u*u/2.0 *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
    return u*u*u*u/4.0 *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
}

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work){

    // compute df/du = df/dR * dR/du

    double direction = dfdu(u)*(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);

    unordered_map_arithmetic(work, du, std::plus<double>(), direction, std::multiplies<double>());

    return 1;
}

double inflow(const vertex& point, const vector<double>& param){

    if (point[1] <= 0.25 || point[1] >= 0.75) {
        return sin(2*M_PI*point[1])*sin(2*M_PI*point[1]);   
    } else {
        return 1;
    }
}

int computeEdgeFlux(vector<double>& edgeflux, double t,
                    const MeshInfo& mi, double ** localvals,
						  const vector<reconstruction>& my_recon,
						  const vector<tensorstencilpoly>& sten_lg,
						  const vector<tensorstencilpoly>& sten_sm){

    int M = mi.MPIglobalCellSize[0];
	 int N = mi.MPIglobalCellSize[1];

    edgeflux.clear();
    edgeflux.resize(M*(N+1) + N*(M+1));
    std::fill(edgeflux.begin(), edgeflux.end(), 0);

    vertexSet constvel {{1.0,0.0},{1.0,0.0},{1.0,0.0}};

    int position = 0;
	 // Need vertical edge only 
    for (int j=0; j<N; j++){
    for (int i=0; i<M-1; i++){

        vertexSet corners = extractCorners(mi, {i,j}); 
        vertexSet edge {corners.at(2), corners.at(1)};

        int neg = j*M+i;
        int pos = j*M+i+1;

        position = M*(N+1) + j*(M+1)+i+1;

        edgeflux.at(position) = advflux_edge(my_recon.at(neg), my_recon.at(pos),
                                sten_lg, sten_sm, localvals, edge, constvel, false,1.0);
    }}

    // inflow boundary with prescribed function
    for (int j=0; j<N; j++){
        vertexSet corners = extractCorners(mi, {0,j});
        vertexSet edge {corners.at(3), corners.at(0)};

        position = M*(N+1) + j*(M+1);

        edgeflux.at(position) = advflux_edge(inflow, {t}, edge, constvel, false, 1.0);

    }

    return 1;
}
