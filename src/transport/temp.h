#ifndef TEMP_H_
#define TEMP_H_

#include "reconstruction.h"
#include "mluse.h"
#include "petsc.h"
//#include "trans_param.h"
#include "util.h"
#include "advectiveflux.h"
#include "diffusiveflux.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

// A temperory file holding simple time stepping for a transport problem
int printSol(int mark, Vec * global, const MeshInfo& mi);

int printGrid(const MeshInfo& mi);

// ====================================================================
double func(const vertex& point, const vector<double>& param);

int RK(double dt, int Nt, Vec * init, const MeshInfo& mi, multilevel& ml, mluse& use, DM dmu, DM dmmesh);

int iRK(double dt, int Nt, Vec * insol, const MeshInfo& mi, multilevel& ml, mluse& use, DM dmu, DM dmmesh, int maxiter);

int getflux(const MeshInfo& mi, multilevel& ml, mluse& use, Vec * now, Vec * flux, DM dmu, DM dmmesh);

int getall(const MeshInfo& mi, multilevel& ml, mluse& use, 
           Vec * now, Vec * flux, Mat *Jacobian, DM dmu, DM dmmesh, 
           const double& dt);

// Output reconstructed values for plotting.
// Compute error at the same time.
double ReconError(const MeshInfo& mi, multilevel& ml, mluse& use,
                  Vec * now, DM dmu, DM dmmesh, const double& h0, int norm);

typedef struct{

    Vec * previous;

    MeshInfo * mi;

    DM dmu;

    DM dmmesh;

    multilevel * ml;

    mluse * use;

    double dt;

    Vec * flux;

} param;

#endif
