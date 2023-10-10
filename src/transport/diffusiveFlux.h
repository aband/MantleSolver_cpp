#ifndef DIFFUSIVEFLUX_H_
#define DIFFUSIVEFLUX_H_

#include "mlwenouse.h"
#include "func.h"

double getDiffusiveFluxInterior(const MLWENOUse& mlu,
                                const MeshInfo& mi,
                                const vertexSet& edge,
                                const vertex& unitNormal,
                                const double& len,
                                const indice& globalCellIn,
                                const indice& globalCellOut,
                                const int& locationIn,
                                const int& locationOut,
                                const valarray<double>& gwe,
                                const valarray<double>& gpe,
                                const double& scale);

#endif
