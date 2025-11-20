#ifndef CLEAN_SOLVER_H_
#define CLEAN_SOLVER_H_

#include <petsc.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include "integral.h"
//#include "input.h"
#include "util.h"

// MFEM parameter header file
#include "myFunc.h"

#include "basis.h"
#include "Hdivmixed.h"
#include "brmixed.h"

// =========================================================================
// The boundary value data structure contains
// 1. global index of degree of freedom and global index of element
// 2. a pair object pairing local degree of freedom and value
struct bndryInfo{
    int        localDOF;
    indice     globalElem;

    double     essenval;
    double     naturval;
};

// Map global indice with boundary values
using bndryVal = std::unordered_map<int, bndryInfo>;

using bndryValGroup = std::unordered_map<std::string, std::vector<bndryVal>>;

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

typedef struct{
    Mat M, Kg, B, Bg, C;
    Vec g, source, neum;
    // Later added vectors
    Vec F, G;
    Vec x, y;
} ReducedSys;

class DarcyStokes{
    public:
        DarcyStokes() {};
        ~DarcyStokes() {};

        int MarkBndryDOFStokes(const MeshInfo& mi,
                               PhysProperty * pp,
                               const std::vector<double>& param);

        int MarkBndryDOFDarcy(const MeshInfo& mi,
                              PhysProperty * pp,
                              const std::vector<double>& param);

    private:

        int AssignLocMatDarcy(LocMat& loc, double theta);
        int AssignLocMatStokes(LocMat& loc, double theta);
        int AssignLocMatCouple(LocMat& loc, double theta);

        

        // Function basis
        basis     * basis_;
        Hdivmixed * hdiv_;
        BRMixed   * br_;

        // Boundary values
        bndryVal bndryStokesEssen_;
        bndryVal bndryStokesNatur_;
        bndryVal bndryDarcyEssen_;
        bndryVal bndryDarcyNatur_;

        /**!
         * Reduced linear system excluding essential boundary conditions
         */
        ReducedSys * reducedDarcy_;
        ReducedSys * reducedStokes_;

        /**!
         * Coupling matrix.
         */
        Mat K;

        /**!
         * Create boundary dof reference mapping
         */
        int bndryDOFStokes_ = 0.0;
        int bndryDOFDarcy_  = 0.0;
 
        int bndryDOFStokesNatur_ = 0.0;
        int bndryDOFDarcyNatur_ = 0.0;

        /**!
         * Reference map tell what kind of boundary condition dof belongs to
         */
        int * refArrayStokesEssen_;
        int * refArrayDarcyEssen_;
 
        unordered_map<int,int> refArrayStokesNatur_;
        unordered_map<int,int> refArrayDarcyNatur_;
 
};

#endif
