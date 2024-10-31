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

vector<double> advFlux(const MeshInfo& mi,
                       const MLWENO::MLWENOUse& mluIn,
                       const MLWENO::MLWENOUse& mluOut,
                       const std::array<vertex,2>& edge,
                       const vertex& unitNormal,
                       const double& len,
                       const indice& globalCellIn,
                       const indice& globalCellOut,
                       const std::string& locationIn,
                       const std::string& locationOut,
                       const vector<double>& LFparam,
                       const vector<double>& direction,
                       const valarray<double>& gpe);

vector<double> advFluxBndry(const MeshInfo& mi,
                            const MLWENO::MLWENOUse& mlu,
                            const std::array<vertex,2>& edge,
                            const vertex& unitNormal,
                            const double& len,
                            const indice& globalCellIn,
                            const int& locationIn,
                            const vector<double>& param,
                            bndryTypeTrans bt);

#endif
