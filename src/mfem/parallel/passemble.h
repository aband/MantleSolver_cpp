#ifndef PASSEMBLE_H_
#define PASSEMBLE_H_

#include <petsc.h>
#include "Hdivmixed.h"
#include "brmixed.h"
#include "util.h"
#include "myFunc.h"
#include "assemble.h"

PetscErrorCode ParallelMatrixAssembleBlock(const MeshInfo& mi, 
                                           basis& basis_,
                                           Hdivmixed& hdiv_,
                                           BRMixed& br_,
                                           PhysProperty * pp,
                                           System * system);

#endif
