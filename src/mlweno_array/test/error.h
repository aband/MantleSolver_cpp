#ifndef ERROR_H_
#define ERROR_H_

#include "util.h"

int printexactsol(const MeshInfo& mi, double t, 
                  double (*func)(const vertex& point,
                                 const vector<double>& param), 
                  int mark, bool grid, const vector<double>& param);

#endif
