#ifndef ADVECTIVEFLUX_H_
#define ADVECTIVEFLUX_H_

#include "mluse.h"
#include "trans_param.h"

// Numerical flux

double edgefluxintegral(const MeshInfo& mi, 
                        const indice& gcellin,
                        const indice& gcellout,
                        const vertexSet& edge,
                        const Tensor<weights>& allwgts,
                        const vector<vertex>& vel,
                        multilevel& ml,
                        mluse& use,
                        double ** lu);

// The most generic function compute edge integral
double edgefluxintegral(const vertexSet& edge,
                        const vector<double>& uin,
                        const vector<double>& uout,
                        const vector<double>& valin,
                        const vector<double>& valout,
                        const vector<double>& dfduin,
                        const vector<double>& dfduout,
                        const vector<vertex>& vel);

// free flow boundary conditions
double edgefluxintegral(const MeshInfo& mi, 
                        const indice& gcell,
                        const vertexSet& edge,
                        const Tensor<weights>& allwgts,
                        const vector<vertex>& vel,
                        multilevel& ml,
                        mluse& use,
                        double ** lu);

// Dirichlet boundary functions
double edgefluxintegral(const vertexSet& edge,
                        const vector<double>& bnval,
                        const vector<vertex>& vel);

double edgefluxintegral(const vertexSet& edge,
                        const double& bnval,
                        const vector<vertex>& vel);

/**!
 * Used in implicit time stepping.
 * Compute Jacobian along with flux
 */
int edgefluxintegral(const MeshInfo& mi, 
                     const indice& gcellin,
                     const indice& gcellout,
                     const vertexSet& edge,
                     const Tensor<weights>& allwgts,
                     const vector<vertex>& vel,
                     multilevel& ml,
                     mluse& use,
                     double ** lu,
                     derivative& der,
                     double& f);

// Used on free outflow boundary condition
int edgefluxintegral(const MeshInfo& mi, 
                     const indice& gcell,
                     const vertexSet& edge,
                     const Tensor<weights>& allwgts,
                     const vector<vertex>& vel,
                     multilevel& ml,
                     mluse& use,
                     double ** lu,
                     derivative& der,
                     double& f);

double getcellflux(const MeshInfo& mi, const indice& gcell,
                   const Tensor<double>& vertedge, 
                   const Tensor<double>& horiedge);

int getcellflux(const MeshInfo& mi, const indice& gcell,
                const Tensor<double>& vertedge, 
                const Tensor<double>& horiedge,
                const Tensor<derivative>& vertedgeder,
                const Tensor<derivative>& horiedgeder,
                double& flux,
                derivative& dflux);

// ======================================================================================

int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts);

int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   Tensor<derivative>& vertedgeder, Tensor<derivative>& horiedgeder,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts);

#endif
