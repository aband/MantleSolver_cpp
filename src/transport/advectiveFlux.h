#ifndef ADVECTIVEFLUX_H_
#define ADVECTIVEFLUX_H_

#include "mlwenouse.h"
#include "func.h"

double getAdvFluxInterior(const MLWENOUse& mlu,
                          const MeshInfo& mi,
                          const vertexSet& edge,
                          const vertex& unitNormal,
                          const double& len,
                          const indice& globalCellIn,
                          const indice& globalCellOut,
                          const int& locationIn,
                          const int& locationout,
                          const valarray<double>& gwe,
                          const valarray<double>& gpe,
                          const double& alpha); 

#endif
