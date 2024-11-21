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

#define COUPLED
#include "passemble.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "bndry.h"
#include "preconst.h"
#include "psolve.h"

// MLWENO parameter header file
#include "mlwenouse.h"
#include "trans_param.h"
#include "advectiveFlux.h"
#include "edgeFlux.h"

extern "C"{
#include "mesh.h"
#include "output.h"
//#include "cgns_io.h"
}

#include<sys/stat.h>

enum transportType {adv, diff, adv_diff, adv_diff_react};

class Driver {
    public:
        //! A constructor
        /**!
         * Construct a driver class.
         * Driver class holding a pointer to meshInfo object.
         * Driver class will be used to interact with underlaying functions.
         */
        Driver() {};

        //! A destructor
        /**!
         * Destruct a Initialize class.
         */
        ~Driver() {delete myPhase;};

        //! A post work clean function
        /**!
         * Clean used dm and vec objects.
         * Should be called at the end of main function.
         */

        int clean();

        /**!
         * MeshInfo struct
         * Can be accessed from outside the class directly.
         */
        MeshInfo mi;

       /**!
         * Data management objects for mesh and solution.
         * showing up in compuation process.
         */ 
        DM dmMesh;       
        DM dmu; 

       /**!
        * Initialize phase package
        */
       Phase * myPhase;

       int CreatePhase();

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

        Vec globalmesh;  

        int PrintMesh();

       /**!
        * Assign Initial cell averaged condition.
        * Specificed for coupled system.
        * global vectors for dimensionless enthalpy and dimensionless Composition.
        */
        Vec globalCD, globalHD;
        Vec localCD, localHD;

        int InitTransport(double (*funcHD)(const valarray<double>& point, const vector<double>& param),
                          double (*funcCD)(const valarray<double>& point, const vector<double>& param));

       // ================================================================================================
       
       /**!
        * Add reconstruction levels to transport problem.
        */
       int AddLevels(const int& stencilSize);
       int AddLevels(const int& stencilSizeX,
                     const int& stencilSizeY);

       int AddLevels(const vector<int>& stencilSizes);

       int AddLevels(const vector<pair<int, int>>& stencilSizes);

       /**!
        * Add default reconstruction levels to transport problem.
        * For advection:
        * (3,2) interior, (3,2,1) on the edge, level 3 on the edge being biased
        * For diffusion:
        * (4,3) interior, (3,2) on the edge, level 3 on the edge being biased
        */
       int PrepareDefaultTransport();

       // ========================================================================

       /**!
        * Create boundary condition vectors
        */
       int PrepareFlow();

       /**!
        * Solve flow at the given time step.
        */
       int SolveFlow(int maxIter, double tolUzawa);

       /**!
        * Scatter distributed vector to all processor.
        * Prepare for velocity reconstruction on gauss points
        */
       int CreateScatterVec();

       /**!
        * Update smoothness indicator and nonlinear weights for all field and mluse
        */
       int UpdateSmoothnessIndicator();

       int UpdateNonlinearWgts();

       /**!
		  * Output of the calculated result
		  */
       int PrintFlow();

       /**!
		  * Print flow during time stepping process
		  */
       int PrintFlowUnscaled(char * filename);

       /**!
        * Print out porosity 
        */
       int PrintPorosity();

       int PrintPorosity(char * filename);

       /**!
        * Print out Pressure
        */
       int PrintPressureConstant();

       int PrintPressureConstantOriginal();

       /**!
        * Show asigned boundary condition.
        * Print assigned boundary conditions type to each dof on the boundary
        */
       int PrintBoundaryDOFs();

       /**!
        * Functions used to compute flux happening on the edges 
        */
       int SingleEdgeFlux(const indice& localedge,
                          extractEdgeInfoFunc edgeinfo, 
                          double& workCD,
                          double& workHD,
                          fluxFunc      fluxfuncAdv, 
                          fluxFuncBndry fluxfuncbndryAdv,
                          fluxFunc      fluxfuncDif,
                          fluxFuncBndry fluxfuncbndryDif);

       int UpdateEdgeFluxAll(vector<double>& edgefluxHD,
                             vector<double>& edgefluxCD,
                             fluxFunc      fluxfuncAdv, 
                             fluxFuncBndry fluxfuncbndryAdv,
                             fluxFunc      fluxfuncDif,
                             fluxFuncBndry fluxfuncbndryDif);

       int ComputeCellFlux(const indice& lCell,
                           double& fluxHD,
                           double& fluxCD,
                           const vector<double>& edgefluxHD,
                           const vector<double>& edgefluxCD);

       // ===================================================================
       /**!
        * counting how many events happen throughout time stepping
        */
       int eventCount;

       char * GetFilename(const char * filename);

    private:

        /**!
         * Record global cell sizes
         */
        int M_;
        int N_;

        /**!
         * Old struct object used in limited functions.
         * Be used for only once.
         */
        MeshParam mp_;

        /**!
         * WENO useage objects
			* Two objects, one for advection and another for diffusion
         */
        MLWENO::MLWENOUse * mluseAdv_;

        MLWENO::MLWENOUse * mluseDif_;

        /**!
         * WENO preparation object.
         */
        MLWENO::MLWENOPrepare * mlpPtr_;

        /**!
			* Position set including all possible positions
			*/
        std::set<std::string> posSet_;
        std::set<std::string> fieldSet_;

        std::unordered_map<std::string, LocFunc> locFuncSet_;

        // ===========================================================

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
        Mat K_;

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

};

// Time stepping struct
// In compleying with C format
typedef struct {

    Driver * driver; 

    double dt;

    int maxIter;

    double tolUzawa;

} ctx_driver;

PetscErrorCode Explicit(TS ts, PetscReal time, Vec U, Vec F, void * ctx);

#endif
