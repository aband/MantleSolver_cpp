#ifndef ERROR_H_
#define ERROR_H_

#include "util.h"
#include "reconstruction.h" 
#include "tensorstencilpoly.h"

#include <petsc.h>
int printexactsol(const MeshInfo& mi, double t, 
                  double (*func)(const vertex& point,
                                 const vector<double>& param), 
                  int mark, bool grid, const vector<double>& param);

int printreconsol(const vector<reconstruction>& my_recon, int M, int N, int mark);

int printreconsol(vector<reconstruction>& my_recon, int M, int N, int mark, 
                  const vector<tensorstencilpoly>& sten_lg,
						const vector<tensorstencilpoly>& sten_sm,
						const MeshInfo& mi,
						double ** locvals);

int printreconsol2(vector<reconstruction>& my_recon, int M, int N, int mark, 
                   const vector<tensorstencilpoly>& sten_lg,
					 	 const vector<tensorstencilpoly>& sten_sm,
					 	 const MeshInfo& mi,
					 	 double ** locvals);

#endif
