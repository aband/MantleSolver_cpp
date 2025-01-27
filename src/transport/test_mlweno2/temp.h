#ifndef TEMP_H_
#define TEMP_H_

#include "reconstruction.h"
#include "petsc.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

// A temperory file holding simple time stepping for a transport problem

double func(const vertex& point, const vector<double>& param);

int RK(double dt, int Nt, Vec * init, const MeshInfo& mi, const multilevel& ml, const mluse& use);

int getflux(const MeshInfo& mi, const multilevel& ml, const mluse& use, Vec * now, Vec * flux);

#endif
