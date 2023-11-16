#ifndef BNDRY_H_
#define BNDRY_H_

#include "myFunc.h"
#include "util.h"
#include "Hdivmixed.h"
#include "brmixed.h"

// The boundary value data structure contains
// 1. global index of degree of freedom and global index of element
// 2. a pair object pairing local degree of freedom and value
using bndryVal = std::unordered_map<int, vector<std::pair<int, double>>>;

typedef struct{
    Vec rhs, source;
    Vec ad, bs, qs, qd;
} RHSVector;

/** !
 * Dirichlet and Neumann boundary conditions are created here.
 *
 * A L2 projection will be used for Dirichlet boundary conditions.
 */

PetscErrorCode AssignValuesRHS(int NS, int ND, int Nelem,
                               Vec * A, Vec * B,
                               RHSVector * rhsv,
                               const bndryVal& bndryStokes,
                               const bndryVal& bndryDarcy);

int MarkBndryDOFStokes(bndryVal& bndryStokes, const MeshInfo& mi, BRMixed& br_);

int MarkBndryDOFDarcy(bndryVal& bndryDarcy, const MeshInfo& mi, Hdivmixed& hdiv_);

int ComputeBndryValsStokes(bndryVal& bndryStokdes, BRMixed& br_,
                           const valarray<double>& gwe,
                           const valarray<double>& gpe);

int ComputeBndryValsDarcy(bndryVal& bndryDarcy, Hdivmixed& hdiv_,
                          const valarray<double>& gwe,
                          const valarray<double>& gpe);

PetscErrorCode CreateRHS(const MeshInfo& mi, 
                         basis& basis_,
                         Hdivmixed& hdiv_,
                         BRMixed& br_,
                         PhysProperty * physproperty,
                         RHSVector * rhsv);
#endif
