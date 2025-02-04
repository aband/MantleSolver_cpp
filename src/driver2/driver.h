#ifndef DRIVER_H_
#define DRIVER_H_

#include <petsc.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include "integral.h"
#include "eutectic.h"
#include "input.h"
#include "util.h"

// MFEM parameter header file
#include "myFunc.h"

//#define COUPLED
#include "passemble.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "bndry.h"
#include "preconst.h"
#include "psolve.h"

// Using new mlweno functions
#include "mluse.h"
#include "reconstruction.h"
#include "stencilpolynomial.h"
#include "tensor.h"

#include "advectiveflux.h"
#include "trans_param.h"

extern "C"{
#include "mesh.h"
#include "output.h"
//#include "cgns_io.h"
}

class Driver {

    public :

        Driver() {};

        ~Driver() { delete myPhase;};

        /**! 
         * Mesh parameters 
         */
        MeshInfo  mi;
        DM dmMesh;
        DM dmu;
        Vec globalmesh;  

       /**!
        * Initialize phase package
        */
       Phase * myPhase;
       int CreatePhase();
       int ShowPhase();
       int withUnit; 

       /**!
        * Create Data management objects.
        * And Mesh vector.
        */
        int CreateMesh(const int& M, const int& N,
                       double L, double H, 
                       double xstart, double ystart,
                       const int& stencilWidthMesh, 
                       const int& stencilWidthU,
                       const bool& physicsScale,
                       const int& meshType); 

        Vec globalCD, globalHD;
        int PrepareTransport(double (*funcHD)(const valarray<double>& point, 
                                              const vector<double>& param),
                             double (*funcCD)(const valarray<double>& point, 
                                              const vector<double>& param));

       /**!
        * Create boundary condition vectors
        */
       int PrepareFlow();

       /**!
        * Solve flow at the given time step.
        */
       int SolveFlow(int maxIter, double tolUzawa, const Tensor<weights>& allwgtsHD, double ** lHD, 
                                                   const Tensor<weights>& allwgtsCD, double ** lCD);

       /**!
        * Scatter distributed vector to all processor.
        * Prepare for velocity reconstruction on gauss points
        */
       int CreateScatterVec();

       /**!
        * Assemble reduced linear system.
        * Reducing boundary dofs
        */
       PetscErrorCode ParallelMatrixAssemble(const Tensor<weights>& allwgtsHD, double ** lHD,
                                             const Tensor<weights>& allwgtsCD, double ** lCD);

       std::vector<double> parameter;

       /**!
        * Interface objects to multileve object
        */
       multilevel ml;

       mluse advection;
       mluse diffusion;

       /**!
        * Simple visualization functions
        */
       int PrintFlowEvent(int mark);

       int PrintPhaseEvent(int mark);

    private:

        /**!
         * Record global cell sizes
         */
        double L_;
        double H_;

        /**!
         * Finite Element spaces.
         */
        basis * basis_;
        Hdivmixed * hdiv_;
        BRMixed * br_;

        /**!
         * Boundary conditions
         */
        // Mark boundary values
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

       /**!
        * Global vector containing results 
        */
       ReducedSys * Result_;

       ScatterResult * sresult_;

       // ==================================================

       int CellAvePorosity(const indice& gcell,
                           const Tensor<weights>& allwgtsHD,
                           double ** lHD,
                           const Tensor<weights>& allwgtsCD,
                           double ** lCD);

       int AssignLocMatStokes(const indice& gcell,
                              const Tensor<weights>& allwgtsHD,
                              double ** lHD,
                              const Tensor<weights>& allwgtsCD,
                              double ** lCD,
                              LocMat * loc);

       int AssignLocMatDarcy(const indice& gcell,
                             const Tensor<weights>& allwgtsHD,
                             double ** lHD,
                             const Tensor<weights>& allwgtsCD,
                             double ** lCD,
                             LocMat * loc);

       int AssignLocMatCouple(const indice& gcell,
                              const Tensor<weights>& allwgtsHD,
                              double ** lHD,
                              const Tensor<weights>& allwgtsCD,
                              double ** lCD,
                              double& k);
};

#endif
