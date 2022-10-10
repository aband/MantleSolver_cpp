#ifndef IMPLICITRK_H_
#define IMPLICITRK_H_

#include<numeric>
#include<assert.h>
#include<algorithm>
#include<petsc.h>

#include "weno_multilevel.h"
#include "flux_multilevel.h"
#include "integral.h"
#include "lapacke.h"

extern "C"{
#include "mesh.h"
}

#include "../test/func.h"


typedef struct{
    DM dm;
    vector<WenoReconstruction *> wr;
    MeshInfo mi;
    int stencil_count;
    int offset;
} Ctx;

PetscErrorCode FormFunction(TS ts, PetscReal time, Vec U, Vec F, void * ctx);

PetscErrorCode FormJacobianIEULER(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void * ctx);

PetscErrorCode SeqImplicitEuler(int stencil_count, vector<double>& linWeights, vector<int *>& rangex, vector<int *>& rangey, MeshInfo& mi, DM dmu, double T, double dt, Vec globalu);


#endif
