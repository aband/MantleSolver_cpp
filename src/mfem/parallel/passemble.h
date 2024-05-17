#ifndef PASSEMBLE_H_
#define PASSEMBLE_H_

#include <petsc.h>
#include "Hdivmixed.h"
#include "brmixed.h"
#include "util.h"
#include "myFunc.h"
#include "assemble.h"
#include "locmat.h"
#include "bndry.h"
#include "solve.h"
#include "shape.h"

PetscErrorCode ParallelAssembleTest();

PetscErrorCode ParallelMatrixAssemble(const MeshInfo& mi,
                                      basis& basis_,
                                      PhysProperty * pp,
                                      const bndryVal& bndryEssenStokes,
                                      ReducedSys * redsysStokes,
                                      const bndryVal& bndryEssenDarcy,
                                      ReducedSys * redsysDarcy,
                                      Mat * K,
                                      BRMixed& br_,
                                      Hdivmixed& hdiv_);

#endif
