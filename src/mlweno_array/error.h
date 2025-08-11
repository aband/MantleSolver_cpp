#ifndef ERROR_H_
#define ERROR_H_

#include "util.h"
#include "reconstruction.h" 
#include "tensorstencilpoly.h"


int printexactsol(const MeshInfo& mi, double t, 
                  double (*func)(const vertex& point,
                                 const vector<double>& param), 
                  int mark, bool grid, const vector<double>& param);

int printreconsol(const vector<reconstruction>& my_recon, int M, int N, int mark);

#endif
