#ifndef FLUX_H_
#define FLUX_H_

#include "weno.h"

solution LaxFriedrichsFlux(const WenoMesh*& wm, int pos, double t,
                           index_set& StencilLarge, vector<index_set>& StencilSmall, point_index& target, 
                           vector<int>& Sorder, vector<int>& Lorder, point target_point);

solution TotalFlux(const WenoMesh*& wm, int pos, double t,
                   index_set& StencilLarge, vector<index_set>& StencilSmall, point_index& target, 
                   vector<int>& Sorder, vector<int>& Lorder,
                   double (*funcX)(valarray<double>& point, const vector<double>& param),
                   double (*funcY)(valarray<double>& point, const vector<double>& param),
                   double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                   double (*dfuncY)(valarray<double>& point, const vector<double>& param));
#endif
