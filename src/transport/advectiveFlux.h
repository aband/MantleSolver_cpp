#ifndef ADVECTIVEFLUX_H_
#define ADVECTIVEFLUX_H_

#include "mlwenouse.h"
#include "trans_param.h"

double getAdvFlux(const MLWENO::MLWENOUse& mlu,
                  const MeshInfo& mi,
                  const std::array<vertex,2>& edge,
                  const vertex& unitNormal,
                  const double& len,
                  const indice& globalCellL,
                  const indice& globalCellR,
                  const int& locationL,
                  const int& locationR,
                  const valarray<double>& gwe,
                  const valarray<double>& gpe,
                  const vector<double>& alpha); 

double getAdvFlux(const MLWENO::MLWENOUse& mlu,
                  const MeshInfo& mi,
                  const std::array<vertex,2>& edge,
                  const vertex& unitNormal,
                  const double& len,
                  const indice& globalCellL,
                  const indice& globalCellR,
                  const int& locationL,
                  const int& locationR,
                  const valarray<double>& gwe,
                  const valarray<double>& gpe,
                  const vector<double>& alpha,
                  bndryTypeTrans bt);

#endif
