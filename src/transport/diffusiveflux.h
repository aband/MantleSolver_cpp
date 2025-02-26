#ifndef DIFFUSIVEFLUX_H_
#define DIFFUSIVEFLUX_H_

#include "mluse.h"
#include "trans_param.h"
#include "lagrange_tmp.h"

double edgefluxintegral(const MeshInfo& mi,
                        const indice& gcellin,
                        const indice& gcellout,
                        const vertexSet& edge,
                        multilevel& ml,
                        mluse& use,
                        double ** lu,
                        const Tensor<weights>& allwgts,
                        const std::string& loc); 
#endif
