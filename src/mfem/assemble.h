#include <petsc.h>
#include "Hdivmixed.h"
#include "brmixed.h"
#include "util.h"
#include "myFunc.h"

typedef struct{
    Mat As,Ad,Bs,Bd,Cs,Cd,K,G,Gp;
    Vec sourceDarcy, sourceStokes;
} System;

typedef struct{
  std::array<double, 64> ad;
  std::array<double, 144> as;
  std::array<double, 8> bd;
  std::array<double, 12> bs;
  std::array<double, 12> rhs;
  double cs;
  double cd;
  double k;

} LocMatrix;

void AssignLocMatrix(const MeshInfo& mi,
                     basis& basis_,
                     Hdivmixed& hdiv_,
                     BRMixed& br_,
                     LocMatrix * locmatrix,
                     PhysProperty * physproperty,
                     const valarray<double>& gwe,
                     const valarray<double>& gpe,
                     const valarray<double>& gwf,
                     const vector<vertex>& gpf);

PetscErrorCode SerialMatrixAssembleBlock(const MeshInfo& mi,
                                         basis& basis_,
                                         Hdivmixed& hdiv_,
                                         BRMixed& br_,
                                         PhysProperty * physpropety,
                                         System * matrix);

PetscErrorCode CreateSchurComplement(System * matrix, 
                                     int nelem, int NS, int ND);
