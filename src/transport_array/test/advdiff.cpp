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

    // Vertical edge flux
    for (int j=0; j<N; j++){
        for (int i=0; i<M-1; i++){

            vertexSet corners = extractCorners(mi, {i,j});
            vertexSet edge    {corners.at(2), corners.at(1)}; 

            int neg = j*M+i;
            int pos = j*M+i+1;

            position = M*(N+1) + j*(M+1)+i+1;
     
            vertexSet quadVel = getQuadVel(velocityField, i, j);;

            edgeflux.at(position) = edgeflux(my_recon.at(neg), my_recon.at(pos),
                                             sten_lg, sten_sm, localvals, edge, );

        }
    }

    return 1;
}
