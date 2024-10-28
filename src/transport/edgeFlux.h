#ifndef EDGEFLUX_H_
#define EDGEFLUX_H_

#include "trans_param.h"

typedef vector<double> (*fluxFunc) (const MeshInfo& mi,
                                    const MLWENO::MLWENOUse& mluIn,
                                    const MLWENO::MLWENOUse& mluOut,
                                    const std::array<vertex,2>& edge,
                                    const vertex& unitNormal,
                                    const double& len,
                                    const indice& globalCellIn,
                                    const indice& globalCellOut,
                                    const std::string& locationIn,
                                    const std::string& locationOut,
                                    const vector<double>& param);

typedef vector<double> (*fluxFuncBndry) (const MeshInfo& mi,
                                         const MLWENO::MLWENOUse& mlu,
                                         const std::array<vertex,2>& edge,
                                         const vertex& unitNormal,
                                         const double& len,
                                         const indice& globalCellIn,
                                         const int& locationIn,
                                         const vector<double>& param,
                                         bndryType bt);

double edgeFlux(const MeshInfo& mi,
                const MLWENO::MLWENOUse& mluIn,
                const MLWENO::MLWENOUse& mluOut,
                const indice& globalCellIn,
                const indice& globalCellOut,
                const std::string& locationIn,
                const std::string& locationOut,
                const std::vector<vertex>& edge,
                const std::vector<vertex>& gauss_p,
                const std::vector<vertex>& velocity,
                const std::valarray<double>& gwe,
                const vetor<double>& param,
                fluxFunc       fluxfunc,
                fluxFuncBndry  fluxfuncbndry);
    
double * edgeFluxAll(const MeshInfo* mi,
                     fluxFunc      fluxfunc,
                     fluxFuncBndry fluxfuncbndry);

#endif
