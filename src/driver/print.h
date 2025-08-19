#ifndef PRINT_H_
#define PRINT_H_

#include "driver.h"

int printCellCenterGrid(const MeshInfo& mi);

int printCellAve(int mrk, Vec * global, const MeshInfo& mi, const char * fieldname);

char * GetFilename(const char * filename, int mark);

int quiverOutputEvent(double * ux, double * uy, double *vx, double * vy, int mark, int M, int N);

int quiverOutputEventTransform(double * ux, double * uy, double *vx, double * vy, int mark, int M, int N, double V0);

int reconPlot(const MeshInfo& mi, multilevel& ml, mluse& use, int mark, Vec * global, bool grid, DM dmu, double h0);

#endif
