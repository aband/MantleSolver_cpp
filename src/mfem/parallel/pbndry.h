#ifndef PBNDRY_H_
#define PBNDRY_H_

#include "myFunc.h"
#include "util.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "bndry.h"

int ParallelMarkBndryDOFs(const MeshInfo& mi,
                          bndryVal& bndryDiri,
                          bndryVal& bndryNeum,
                          basis& basis_,
                          PhysProperty * pp,
                          BRMixed& br_);

#endif
