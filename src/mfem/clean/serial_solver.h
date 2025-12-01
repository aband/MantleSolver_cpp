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

        int init(const MeshInfo& mi,
                 PhysProperty * pp,
                 const std::vector<double>& param);

        int printBndryAll();

    private:

        int ComputeEssenBndryAll(const MeshInfo& mi,
                                 PhysProperty * pp,
                                 const std::vector<double>& param);

        int ComputeNaturBndryAll(const MeshInfo& mi,
                                 PhysProperty * pp,
                                 const std::vector<double>& param);

        double AssignBndrySupVal(const vertexSet& edgeCorner,
                                 const vertex& nu,
                                 PhysProperty * pp);

        std::array<double, 2> AssignBndryValsDarcy(const vertexSet& edgeCorner,
                                                   const vertex& nu,
                                                   const double& len,
                                                   const int& edge,
                                                   PhysProperty * pp);

        int computeEssenVals(const MeshInfo& mi, int i, int j, int edge, PhysProperty * pp);

        int AssignLocMatDarcy(LocMat& loc, double theta);
        int AssignLocMatStokes(LocMat& loc, double theta);
        int AssignLocMatCouple(LocMat& loc, double theta);

        // normal dof 
        std::set<int> top_normal {11,7,6};
        std::set<int> left_normal {0,3,8};
        std::set<int> right_normal {1,2,10};
        std::set<int> bottom_normal {4,5,9};

        // tangent dof
        std::set<int> top_tang {3,2};
        std::set<int> left_tang {7,4};
        std::set<int> right_tang {5,6};
        std::set<int> bottom_tang {0,1};

        // Function basis
        basis     basis_;
        Hdivmixed hdiv_;
        BRMixed   br_;

        // Cleared boundary values 
        bndryVal bndryStokesEssen_;
        bndryVal bndryStokesNatur_;
        bndryVal bndryDarcyEssen_;
        bndryVal bndryDarcyNatur_;

        /* ====================================================  
         Boundary values containing four edges numberred 0-3
         The numbering order is the same as the element edge number
		   0 -- left
		   1 -- bottom
		   2 -- right
		   3 -- top
			The corner dofs are defined and dealt with separately
		  ==================================================== */ 

        bndryVal bndryStokesEssenAll;
        bndryVal bndryStokesNaturAll;
        bndryVal bndryDarcyEssenAll;
        bndryVal bndryDarcyNaturAll;

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
