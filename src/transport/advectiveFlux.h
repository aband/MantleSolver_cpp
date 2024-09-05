#ifndef ADVECTIVEFLUX_H_
#define ADVECTIVEFLUX_H_

#include "mlwenouse.h"
#include "trans_param.h"

// Denoting different boundary types for a given physical domain
enum bndryType {"wall", "free", "dirichlet", "neumann", "absorb", "reflect", "periodic", "flux"};

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
                  const double& alpha); 

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
                  const double& alpha,
                  bndryType bt);

#endif
