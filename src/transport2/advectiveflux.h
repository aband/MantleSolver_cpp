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

int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts);

double getcellflux(const MeshInfo& mi, const indice& gcell,
                   const Tensor<double>& vertedge, 
                   const Tensor<double>& horiedge);

#endif
