#ifndef DIFFUSIVEFLUX_H_
#define DIFFUSIVEFLUX_H_

#include "mlwenouse.h"
#include "trans_param.h"

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

double diffFlux(const valarray<double>& gwe,
                const vector<vertex>& param,
                const vector<double>& uIn,
                const vector<double>& uOut,
                const vertex& unitnormal,
                const double& len);

double diffFluxBndry(const MeshInfo& mi,
                     const valarray<double>& gwe,
                     const vector<vertex>& vel,
                     const vector<double>& u,
                     const vertex& unitnormal,
                     const double& len,
                     const indice& gCell,
                     const int& edgeflag,
                     const std::string& field);

#endif
