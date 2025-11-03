#ifndef DIFFUSIVEFLUX_ARRAY_H_
#define DIFFUSIVEFLUX_ARRAY_H_

#include "tensorstencilpoly.h"
#include "reconstruction.h"

#include "lagrange_tmp.h"

#include "transfunc.h"

double diffflux_edge(const reconstruction& recon_neg,
                     const reconstruction& recon_pos,
                     const vector<tensorstencilpoly>& sten_lg,
                     const vector<tensorstencilpoly>& sten_sm,
                     double ** localvals,
                     const vertexSet& edge);

#endif
