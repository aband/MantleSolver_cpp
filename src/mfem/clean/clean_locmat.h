#ifndef CLEAN_LOCMAT_H_
#define CLEAN_LOCMAT_H_

#include <petsc.h>
#include "Hdivmixed.h"
#include "brmixed.h"
#include "util.h"
#include "shape.h"
#include "myFunc.h"

typedef struct{

  std::vector<double> A;
  std::vector<double> B;
  double C;
  std::vector<double> f;

} LocMat;

typedef struct{

  std::vector<double> edgeporo;
  std::vector<double> cellporo;
  double aveporo;

} poroSet;

int AssignLocMatStokes(const MeshInfo& mi,
                       BRMixed& br_,
                       basis& basis_,
                       LocMat& loc,
                       double theta, 
                       const poroSet& poro);

int AssignLocMatDarcy(const MeshInfo& mi,
                      Hdivmixed& hdiv_,
                      basis& basis_,
                      LocMat& loc,
                      double theta,
                      const poroSet& poro);

int AssignLocMatCouple(const MeshInfo& mi,
                       );

#endif
