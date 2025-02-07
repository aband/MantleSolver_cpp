#ifndef EDGEFLUX_H_
#define EDGEFLUX_H_

#include "trans_param.h"
#include "advectiveFlux.h"

typedef double (*fluxFunc) (const valarray<double>& gwe,
                            const vector<vertex>& vel,
                            const vector<double>& uIn,
                            const vector<double>& uOut,
                            const vertex& unitnormal,
                            const double& len);

typedef double (*fluxFuncBndry) (const MeshInfo& mi,
                                 const valarray<double>& gwe,
                                 const vector<vertex>& vel,
                                 const vector<double>& u,
                                 const vertex& unitnormal,
                                 const double& len,
                                 const indice& gCell,
                                 const int& edgeflag,
                                 const std::string& field);

#endif
