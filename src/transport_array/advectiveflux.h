#ifndef ADVECTIVEFLUX_ARRAY_H_
#define ADVECTIVEFLUX_ARRAY_H_

#include "tensorstencilpoly.h"
#include "reconstruction.h"

#include "transfunc.h"

inline double LFflux(double fneg, double fpos, double uneg, double upos, double alpha){
    return 0.5*(fneg + fpos - alpha*(upos-uneg));    
}

double advflux_edge(const reconstruction& recon_neg,
					     const reconstruction& recon_pos,
                    const vector<tensorstencilpoly>& sten_lg,
                    const vector<tensorstencilpoly>& sten_sm,
                    double ** localvals,
						  const vertexSet& edge,
						  const vector<vertex>& vel,
						  bool localLF,
						  double gLF);

double advflux_edge(const reconstruction& recon,
                    const vector<tensorstencilpoly>& sten_lg,
                    const vector<tensorstencilpoly>& sten_sm,
                    double ** localvals,
						  const vertexSet& edge,
						  const vector<vertex>& vel,
						  bool localLF,
						  double gLF);

double advflux_edge(double (*func)(const vertex& point,
								           const vector<double>& param),
					     const vector<double>& param,
					     const vertexSet& edge,
						  const vector<vertex>& vel,
						  bool localLF,
						  double gLF);

int testlink();

#endif
