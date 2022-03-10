#ifndef FLUX_H_
#define FLUX_H_

#include "weno.h"

solution LaxFriedrichsFlux(const WenoMesh*& wm, int pos, double t, WenoReconst**& wr, 
                           point target_point, point_index& target, 
                           double (*funcX)(valarray<double>& point, const vector<double>& param),
                           double (*funcY)(valarray<double>& point, const vector<double>& param),
                           double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                           double (*dfuncY)(valarray<double>& point, const vector<double>& param));

solution TotalFlux(const WenoMesh*& wm, int pos, double t,
                   point_index& target, WenoReconst**& wr,
                   double (*funcX)(valarray<double>& point, const vector<double>& param),
                   double (*funcY)(valarray<double>& point, const vector<double>& param),
                   double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                   double (*dfuncY)(valarray<double>& point, const vector<double>& param));


// The following functions recalculate weno basis polynomials each iteration
// Do not use these functions directly.
// Use as benchmark instead.
solution LaxFriedrichsFlux(const WenoMesh*& wm, int pos, double t,
                           index_set& StencilLarge, vector<index_set>& StencilSmall, point_index& target, 
                           vector<int>& Sorder, vector<int>& Lorder, point target_point, 
                           double (*funcX)(valarray<double>& point, const vector<double>& param),
                           double (*funcY)(valarray<double>& point, const vector<double>& param),
                           double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                           double (*dfuncY)(valarray<double>& point, const vector<double>& param));

solution TotalFlux(const WenoMesh*& wm, int pos, double t,
                   index_set& StencilLarge, vector<index_set>& StencilSmall, point_index& target, 
                   vector<int>& Sorder, vector<int>& Lorder,
                   double (*funcX)(valarray<double>& point, const vector<double>& param),
                   double (*funcY)(valarray<double>& point, const vector<double>& param),
                   double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                   double (*dfuncY)(valarray<double>& point, const vector<double>& param));
#endif
