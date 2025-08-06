#ifndef ADVECTIVEFLUX_H_
#define ADVECTIVEFLUX_H_

#include "tensorstencilpoly.h"
#include "reconstruction.h"

inline double LFflux(double fneg, double fpos, double uneg, double upos, double alpha){
    return 0.5*(fneg + fpos - alpha*(upos-uneg));    
}

double edgefluxintegral(const reconstruction& recon_out,
                        const reconstruction& recon_in);


#endif
