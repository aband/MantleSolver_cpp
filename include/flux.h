#ifndef FLUX_H_
#define FLUX_H_

#include "weno.h"

solution LaxFriedrichsFlux(const WenoMesh*& wm, solution u_in, int pos, double t,
                           index_set& StencilLarge, vector<index_set>& StencilSmall, point_index& target, 
                           vector<int>& Sorder, vector<int>& Lorder, point target_point);

#endif
