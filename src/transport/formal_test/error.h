#ifndef MLWENOERROR_H_
#define MLWENOERROR_H_

#include "mluse.h"

#include "util.h"
#include "advectiveflux.h"
#include "diffusiveflux.h"

int reconPlot(const MeshInfo& mi, multilevel& ml, mluse& use, int mark, Vec * global, bool grid, DM dmu, double h0);

int simpleRK(double dt, int Nt, Vec * insol, const MeshInfo& mi, multilevel& ml,  mluse& use, DM dmu, DM dmmesh);

int simpleSSP2RK(double dt, int Nt, Vec * insol, const MeshInfo& mi, multilevel& ml,  mluse& use, DM dmu, DM dmmesh);

int simpleSSP3RK(double dt, int Nt, Vec * insol, const MeshInfo& mi, multilevel& ml,  mluse& use, DM dmu, DM dmmesh);

// Compute effective reconstruction order
int eff_order(const MeshInfo& mi, multilevel& ml, mluse& use, int mark, Vec * global, bool grid, DM dmu, double h0);

#endif
