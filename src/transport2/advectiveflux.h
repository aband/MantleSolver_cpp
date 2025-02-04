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
                        multilevel& ml,
                        mluse& use,
                        double ** lu);

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
