#ifndef BNDRY_H_
#define BNDRY_H_

#include "myFunc.h"
#include "util.h"
#include "Hdivmixed.h"
#include "brmixed.h"

// The boundary value data structure contains
// 1. global index of degree of freedom and global index of element
// 2. a pair object pairing local degree of freedom and value
struct bndryInfo{
    int    localDOF;
    double DirichletVal;  
    indice globalElem;
};

using bndryVal = std::unordered_map<int, bndryInfo>;

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

/** !
 * Functions mark Dirichlet boundary values.
 * Darcy and Stokes boundary values are assigned differently.
 */
int MarkBndryDOFStokes(bndryVal& bndryStokes, 
                       const MeshInfo& mi, 
                       basis& basis_,
                       BRMixed& br_);

int MarkBndryDOFDarcy(bndryVal& bndryDarcy, 
                      const MeshInfo& mi,
                      basis& basis_,
                      Hdivmixed& hdiv_);

/** !
 * Compute Dirichlet values locally.
 * 1. Dirichlet values are assigned to Stokes part directly.
 * 2. Dirichlet values are assigned to Darcy part via L2 projection.
 */
std::array<double, 2> AssignBndryValsDarcy(const indice& global, 
                                           const int& edge, 
                                           basis& basis_,
                                           Hdivmixed& hdiv_,
                                           const vertexSet& edgeCorner,
                                           const double& len,
                                           const valarray<double>& gwe,
                                           const valarray<double>& gpe);

double AssignBndrySupVal(const vertexSet& edgeCorner, 
                         const vertex& nu,
                         const valarray<double>& gwe,
                         const valarray<double>& gpe);

/** !
 * Assign Neumann or Dirichlet boundary to different elements.
 * A Neumann Dirichlet mixed boundary condition.
 * Only two different kinds of boundary conditions.
 */
bool Is_Dirichlet(const indice& global);

/** !
 * Create Matrix Kg and Vector g regarding Dirichlet boundary condition
 * both serial and parallel versions of functions are provided.
 */
petscErrorCode CreateDirichletMatVecSerial();


PetscErrorCode CreateRHS(const MeshInfo& mi, 
                         basis& basis_,
                         Hdivmixed& hdiv_,
                         BRMixed& br_,
                         PhysProperty * physproperty,
                         RHSVector * rhsv);
#endif
