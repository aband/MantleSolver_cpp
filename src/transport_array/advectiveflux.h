#ifndef ADVECTIVEFLUX_H_
#define ADVECTIVEFLUX_H_

#include "tensorstencilpoly.h"
#include "reconstruction.h"

inline double LFflux(double fneg, double fpos, double uneg, double upos, double alpha){
    return 0.5*(fneg + fpos - alpha*(upos-uneg));    
}

int advflux_all(const vector<reconstruction>& my_recon,
                const vector<tensorstencilpoly>& sten_lg,
                const vector<tensorstencilpoly>& sten_sm,
                double ** localvals,
                vector<double>& advflux,
                const MeshInfo& mi);

#endif
