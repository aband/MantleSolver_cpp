#ifndef ADVECTIVEFLUX_H_
#define ADVECTIVEFLUX_H_

#include "mlwenouse.h"
#include "trans_param.h"

vector<double> advFlux(const MeshInfo& mi,
                       const MLWENO::MLWENOUse& mlu,
                       const vector<vertex>& edge,
                       const vertex& unitNormal,
                       const double& len,
                       const indice& gCellIn,
                       const indice& gCellOut,
                       const std::string& locIn,
                       const std::string& locOut,
                       const vector<double>& LFparam,
                       const vector<double>& direction,
                       const valarray<double>& gpe);

vector<double> advFluxBndry(const MeshInfo& mi,
                            const MLWENO::MLWENOUse& mlu,
                            const vector<vertex>& edge,
                            const vertex& unitNormal,
                            const double& len,
                            const indice& gCell,
                            const std::string& loc,
                            const vector<double>& direction,
                            const vector<double>& LFparam,
                            const valarray<double>& gpe,
                            const bndryTypeTrans& bt,
                            const int& flag,
                            const int& locedge);

double advFlux(const valarray<double>& gwe,
               const vector<vertex>& velOut, 
               const vector<vertex>& velIn,
               const vector<double>& uIn, 
               const vector<double>& uOut,
               const vertex& unitnormal,
               const double& len);

double advFluxBndry(const valarray<double>& gwe,
                    const vector<vertex>& vel,
                    const vector<double>& u,
                    const vertex& unitnormal,
                    const double& len,
                    const bndryTypeTrans& bt,
                    const int& edgetype);

#endif
