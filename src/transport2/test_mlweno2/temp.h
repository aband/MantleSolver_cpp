#ifndef TEMP_H_
#define TEMP_H_

#include "reconstruction.h"
#include "mluse.h"
#include "petsc.h"
#include "trans_param.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

// A temperory file holding simple time stepping for a transport problem

double func(const vertex& point, const vector<double>& param);

int RK(double dt, int Nt, Vec * init, const MeshInfo& mi, multilevel& ml, const mluse& use, DM dmu, DM dmmesh);

int getflux(const MeshInfo& mi, multilevel& ml, mluse& use, Vec * now, Vec * flux, DM dmu, DM dmmesh);

#endif
