#ifndef FLUX_MULTILEVEL_H_
#define FLUX_MULTILEVEL_H_

#include "weno_multilevel.h"

double LaxFriedrichsFlux(MeshInfo* mi, int pos, double t, WenoReconstruction**& wr,
                         point& target_point, point_index& target_index,
                         double (*funcX)(valarray<double>& point, const vector<double>& param),
                         double (*funcY)(valarray<double>& point, const vector<double>& param),
                         double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                         double (*dfuncY)(valarray<double>& point, const vector<double>& param));

double TotalFlux(MeshInfo* mi, int pos, double t,
                 point_index& target_index, WenoReconstruction**& wr,
                 double (*funcX)(valarray<double>& point, const vector<double>& param),
                 double (*funcY)(valarray<double>& point, const vector<double>& param),
                 double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                 double (*dfuncY)(valarray<double>& point, const vector<double>& param));

#endif
