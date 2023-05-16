#ifndef TIMESTEPPING_H_
#define TIMESTEPPING_H_

#include "transport.h"

/**
 * Define struct object Ctx.
 */
typedef struct {
    transport * trPtr;
    DM dmu;
    MeshInfo * mi;
} Ctx;

/**
 * Explicit time stepping.
 * Using time stepping object provided by Petsc.
 */
PetscErrorCode Explicit(TS ts, PetscReal time, Vec U, Vec F, void* ctx);

PetscErrorCode ExplicitAdvection(TS ts, PetscReal time, Vec U, Vec F, void* ctx);

PetscErrorCode ExplicitDiffusion(TS ts, PetscReal time, Vec U, Vec F, void* ctx);

/**
 * Implicit time stepping.
 * Using time stepping object provided by Petsc.
 */
PetscErrorCode FormJacobian(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void* ctx);

PetscErrorCode FormJacobianAdvection(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void* ctx);

PetscErrorCode FormJacobianDiffusion(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void* ctx);

#endif
