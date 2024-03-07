#ifndef BNDRY_H_
#define BNDRY_H_

#include "myFunc.h"
#include "util.h"
#include "Hdivmixed.h"
#include "brmixed.h"

enum bndryType {dirichlet, neumann, robin};

// The boundary value data structure contains
// 1. global index of degree of freedom and global index of element
// 2. a pair object pairing local degree of freedom and value
struct bndryInfo{
    int        localDOF;
    double     DirichletVal;  
    indice     globalElem;
    bndryType  bt;
};

using bndryVal = std::unordered_map<int, bndryInfo>;

typedef struct{
    Mat M, Kg, B, Bg;
    Vec g, source;
} ReducedSys;

/** !
 * Dirichlet and Neumann boundary conditions are created here.
 *
 * A L2 projection will be used for Dirichlet boundary conditions.
 */

/** !
 * Functions mark Dirichlet boundary values.
 * Darcy and Stokes boundary values are assigned differently.
 */
int MarkBndryDOFStokes(bndryVal& bndryStokes, 
                       const MeshInfo& mi, 
                       basis& basis_,
                       BRMixed& br_,
                       PhysProperty * pp);

int MarkBndryDOFDarcy(bndryVal& bndryDarcy, 
                      const MeshInfo& mi,
                      basis& basis_,
                      Hdivmixed& hdiv_,
                      PhysProperty * pp);

/** !
 * Compute Dirichlet values locally.
 * 1. Dirichlet values are assigned to Stokes part directly.
 * 2. Dirichlet values are assigned to Darcy part via L2 projection.
 */
std::array<double, 2> AssignBndryValsDarcy(const indice& global, 
                                           const int& edge, 
                                           basis& basis_,
                                           Hdivmixed& hdiv_,
                                           PhysProperty * pp,
                                           const vertexSet& edgeCorner,
                                           const double& len,
                                           const valarray<double>& gwe,
                                           const valarray<double>& gpe);

double AssignBndrySupVal(const vertexSet& edgeCorner, 
                         const vertex& nu,
                         const valarray<double>& gwe,
                         const valarray<double>& gpe,
                         PhysProperty * pp);

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
PetscErrorCode CreateReducedSerial(ReducedSys * reducedsys,
                                   Mat * fullM, Mat * fullB,
                                   Vec * fullSource,
                                   const bndryVal& bndryval);

#endif
