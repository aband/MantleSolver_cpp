// A advection diffusion test case

#include "transfunc.h"

// Dirichlet inflow boundary condition
double inflow(const vertex& point, const vector<double>& param){

    return 0.0;
}

int computeEdgeFlux(const vector<vertex>& velocityField,
                    vector<double>& edgeflux, double t,
                    const MeshInfo& mi, double ** localvals,
                    const vector<reconstruction>& my_recon,
                    const vector<tensorstencilpoly>& sten_lg,
                    const vector<tensorstencilpoly>& sten_sm){

    edgeflux.clear();
    edgeflux.resize(M*(N+1) + N*(M+1));

    std::fill(edgeflux.begin(), edgeflux.end(), 0);


    return 1;
}
