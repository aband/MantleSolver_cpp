#ifndef EDGEFLUX_H_
#define EDGEFLUX_H_

#include "trans_param.h"
#include "advectiveFlux.h"

typedef vector<double> (*fluxFunc) (const MeshInfo& mi,
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

typedef vector<double> (*fluxFuncBndry) (const MeshInfo& mi,
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

template <typename T>
struct edgeEnds{
    T start;
    T end;
};

double * edgeFluxAll(const MeshInfo* mi,
                     fluxFunc      fluxfunc,
                     fluxFuncBndry fluxfuncbndry);

#endif
