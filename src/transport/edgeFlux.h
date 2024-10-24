#ifndef EDGEFLUX_H_
#define EDGEFLUX_H_

typedef double (*fluxFunc) (const MLWENO::MLWENOUse& mlu,
                            const MLWENO::MLWENOUse& mlu,
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
                            const vector<double>& alpha);

typedef double (*fluxFuncBndry) (const MLWENO::MLWENOUse& mlu,
                                 const MeshInfo& mi,
                                 const std::array<vertex,2>& edge,
                                 const vertex& unitNormal,
                                 const double& len,
                                 const indice& globalCellIn,
                                 const int& locationIn,
                                 const valarray<double>& gwe,
                                 const valarray<double>& gpe,
                                 const vector<double>& alpha,
                                 bndryTypeAdv bt);

double * edgeFluxAll(const MeshInfo* mi,
                     fluxFunc      fluxfunc,
                     fluxFuncBndry fluxfuncbndry);

#endif
