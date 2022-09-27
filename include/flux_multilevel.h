#ifndef FLUX_MULTILEVEL_H_
#define FLUX_MULTILEVEL_H_

#include "weno_multilevel.h"

double LaxFriedrichsFlux(const MeshInfo& mi, int pos, double t, vector<WenoReconstruction*>& wr,
                         point& target_point, point_index& target_index,
                         double (*funcX)(valarray<double>& point, const vector<double>& param),
                         double (*funcY)(valarray<double>& point, const vector<double>& param),
                         double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                         double (*dfuncY)(valarray<double>& point, const vector<double>& param));

double TotalFlux(const MeshInfo& mi, int pos, double t,
                 point_index& target_index, vector<WenoReconstruction*>& wr,
                 double (*funcX)(valarray<double>& point, const vector<double>& param),
                 double (*funcY)(valarray<double>& point, const vector<double>& param),
                 double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                 double (*dfuncY)(valarray<double>& point, const vector<double>& param));

vector<double> DerivLaxFriedrichFlux(const MeshInfo& mi, double t,
                                     point_index& target_index, vector<WenoReconstruction*>& wr,
                                     double (*funcX)(valarray<double>& point, const vector<double>& param),
                                     double (*funcY)(valarray<double>& point, const vector<double>& param),
                                     double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                                     double (*dfuncY)(valarray<double>& point, const vector<double>& param));

#endif
