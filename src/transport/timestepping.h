#ifndef TIMESTEPPING_H_
#define TIMESTEPPING_H_

#include "transport.h"

/**
 * Define struct object Ctx.
 */
typedef struct {
    transport* trPtr;
    DM dmu;
    MeshInfo * mi;
} Ctx;

/**
 * Explicit Eurler.
 * Using time stepping object provided by Petsc.
 */
PetscErrorCode Explicit(TS ts, PetscReal time, Vec U, Vec F, void* ctx);

#endif
