#ifndef ADVECTIVEFLUX_H_
#define ADVECTIVEFLUX_H_

#include "mlwenouse.h"
#include "trans_param.h"

double advFlux(const valarray<double>& gwe,
               const vector<vertex>& vel, 
               const vector<double>& uIn, 
               const vector<double>& uOut,
               const vertex& unitnormal,
               const double& len);

double advFluxBndry(const MeshInfo& mi,
                    const valarray<double>& gwe,
                    const vector<vertex>& vel,
                    const vector<double>& u,
                    const vertex& unitnormal,
                    const double& len,
                    const indice& gCell,
                    const int& edgeflag,
                    const std::string& field);

#endif
