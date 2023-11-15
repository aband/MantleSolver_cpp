#ifndef BNDRY_H_
#define BNDRY_H_

// The boundary value data structure contains
// 1. global index of degree of freedom
// 2. a pair object pairing local degree of freedom and value
using bndryVal = std::unordered_map<int, std::pari<int, double>>;

/** !
 * Dirichlet and Neumann boundary conditions are created here.
 *
 * A L2 projection will be used for Dirichlet boundary conditions.
 */

double


#endif
