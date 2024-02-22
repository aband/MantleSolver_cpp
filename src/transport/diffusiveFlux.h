#ifndef DIFFUSIVEFLUX_H_
#define DIFFUSIVEFLUX_H_

#include "mlwenouse.h"
#include "../driver/mlwenoTest/transportTest/myfunc.h"

double getDifFluxInterior(const MLWENO::MLWENOUse& mlu,
                          const MeshInfo& mi,
                          const std::array<vertex,2>& edge,
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
