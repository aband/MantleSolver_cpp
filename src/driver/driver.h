#ifndef DRIVER_H_
#define DRIVER_H_

#include <petsc.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include "integral.h"
//#include "eutectic.h"
#include "eutectic_rescaled.h"
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

        int PrepareTransport(const std::vector<double>& restartHD,
                             const std::vector<double>& restartCD,
                             double (*funcHD)(const valarray<double>& point, 
                                              const vector<double>& param),
                             double (*funcCD)(const valarray<double>& point, 
                                              const vector<double>& param));

        double h0;
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

       double HDbottom;
       double CDbottom;

        /**!
         * Time stepping function
         */
       int RK(double dt, double Tmax, int maxIter, double tolUzawa);

       int RK_Pause(double dt, double Tmax, int maxIter, 
                    double tolUzawa, int interval);

       int getFluxAll(Vec * fCD, Vec * fHD,  Vec * gCD, Vec * gHD,
                      int t, double dt, int maxIter, double tolUzawa, int interval);

       int getFluxAll(const Tensor<vertexSet>& phasevel_vert, 
                      const Tensor<vertexSet>& phasevel_hori, 
                      const Tensor<vertexSet>& effvel_vert, 
                      const Tensor<vertexSet>& effvel_hori, 
                      const Tensor<vertexSet>& solidvel_vert, 
                      const Tensor<vertexSet>& solidvel_hori,
                      Vec * fCD, Vec * fHD, Vec * gCD, Vec * gHD, 
                      int t, double dt, int maxIter, double tolUzawa, int interval);

       int SSP2RK(double dt, double Tmax, int maxIter, 
                  double tolUzawa, int interval);

       int SSP2RK_Pause(double dt, double Tmax, int maxIter, 
                        double tolUzawa, int interval);

       /**!
        * Simple visualization functions
        */
       int PrintFlowEvent(int mark);
       int PrintFlowEventTransform(int mark);

       int PrintPhaseEvent(int mark);

       int PrintPressureSerialApprox(int mark);

       int PrintEffVel(int mark, int side,
                       const Tensor<weights>& allwgtsHD, double ** lHD,
                       const Tensor<weights>& allwgtsCD, double ** lCD);

       int PrintTensorVel(const Tensor<vertexSet>& pointvel,
                          int mark,
                          const char * filedname);

       double V0;

       int start;

       // Case study ===============================================================
       int computeEffVel_case(const vector<vertex>& gaussp,
                              const vertexSet& edgep,
                              const indice& gcellin, const indice& gcellout,
                              const Tensor<weights>& allwgts, double ** lphi,
                              vector<vertex>& vel);
 
       int computeEffVel_case(const vector<vertex>& gaussp,
                              const vertexSet& edgep,
                              const indice& gcell,
                              const Tensor<weights>& allwgts, double ** lphi,
                              vector<vertex>& vel);

       int updateEdgeFlux_case(Tensor<double>& vertedge, Tensor<double>& horiedge,
                                const Tensor<weights>& allwgts, double ** lphi);

       int CellAvePorosity_case(const indice& gcell, 
                                const Tensor<weights>& allwgts,
                                double ** lphi);

       int AssignLocMatStokes_case(const indice& gcell,
                                   const Tensor<weights>& allwgts,
                                   double ** lphi,
                                   LocMat * loc);

       int AssignLocMatDarcy_case(const indice& gcell,
                                  const Tensor<weights>& allwgts,
                                  double ** lphi,
                                  LocMat * loc);
 
       int AssignLocMatCouple_case(const indice& gcell,
                                   const Tensor<weights>& allwgts,
                                   double ** lphi,
                                   double& k);

       int PrepareTransport_case(double (*func)(const valarray<double>& point,
                                                const vector<double>& param) );

       int ParallelMatrixAssemble_case(const Tensor<weights>& allwgts,
                                       double ** lphi);

       int RK_case(double dt, double Tmax, int maxIter, 
                   double tolUzawa);

       int SolveFlow_case(int maxIter, double tolUzawa, 
                          const Tensor<weights>& allwgts, double ** lphi);

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

       // =======================================================
       int computephase(const std::vector<vertex>& gaussp,
                        const vertexSet& edgep,
                        const indice& gcell,
                        const Tensor<weights>& allwgtsHD, double ** lHD,
                        const Tensor<weights>& allwgtsCD, double ** lCD,
                        vector<double>& cs,
                        vector<double>& cl,
                        vector<double>& phif,
                        vector<double>& TD,
                        vector<double>& dTdH,
								vector<double>& CD,
								vector<double>& HD);

       // Function used in the interior region
       int computeEffVel(const vector<vertex>& gaussp,
                         const vertexSet& edgep,
                         const indice& gcellin, const indice& gcellout,
                         const Tensor<weights>& allwgtsHD, double ** lHD,
                         const Tensor<weights>& allwgtsCD, double ** lCD,
                         vector<vertex>& effvelHD, 
                         vector<vertex>& effvelCD);

       int computeEffVel(const vector<vertex>& gaussp,
                         const vertexSet& edgep,
                         const indice& gcellin, const indice& gcellout,
                         const Tensor<weights>& allwgtsHD, double ** lHD,
                         const Tensor<weights>& allwgtsCD, double ** lCD,
                         vector<vertex>& effvel,
                         vector<vertex>& phasevel,
                         vector<vertex>& solidvel,
                         vector<double>& TDin,
                         vector<double>& TDout,
                         vector<double>& dTdHin,
                         vector<double>& dTdHout,
								 vector<double>& CDin,
								 vector<double>& CDout,
								 vector<double>& HDin,
								 vector<double>& HDout);

       // Function used on the boundary
       int computeEffVel(const vector<vertex>& gaussp,
                         const vertexSet& edgep,
                         const indice& gcell,
                         const Tensor<weights>& allwgtsHD, double ** lHD,
                         const Tensor<weights>& allwgtsCD, double ** lCD,
                         vector<vertex>& effvelHD, 
                         vector<vertex>& effvelCD);

       int computeEffVel(const vector<vertex>& gaussp,
                         const vertexSet& edgep,
                         const indice& gcell,
                         const Tensor<weights>& allwgtsHD, double ** lHD,
                         const Tensor<weights>& allwgtsCD, double ** lCD,
                         vector<vertex>& effvel,
                         vector<vertex>& phasevel,
                         vector<vertex>& solidvel,
                         vector<double>& TDin,
                         vector<double>& dTdHin,
								 vector<double>& CD,
								 vector<double>& HD);

       int computeFaceVel(const vector<vertex>& gaussp,
                          const indice& gcell, 
                          const Tensor<weights>& allwgtsHD, double ** lHD,
                          const Tensor<weights>& allwgtsCD, double ** lCD,
                          vector<vertex>& phasevel,
                          vector<double>& TD);

       int updateEdgeFlux(Tensor<double>& vertedgeHD, Tensor<double>& horiedgeHD,
                          Tensor<double>& vertedgeCD, Tensor<double>& horiedgeCD,
                          const Tensor<weights>& allwgtsHD, double ** lHD,
                          const Tensor<weights>& allwgtsCD, double ** lCD);

       int updateEdgeFlux(const Tensor<vertexSet>& phasevel_vert, 
                          const Tensor<vertexSet>& phasevel_hori, 
                          const Tensor<vertexSet>& effvel_vert, 
                          const Tensor<vertexSet>& effvel_hori, 
                          const Tensor<vertexSet>& solidvel_vert, 
                          const Tensor<vertexSet>& solidvel_hori,
                          Tensor<double>& vertedgeHD, Tensor<double>& horiedgeHD,
                          Tensor<double>& vertedgeCD, Tensor<double>& horiedgeCD,
                          const Tensor<weights>& allwgtsHD, double ** lHD,
                          const Tensor<weights>& allwgtsCD, double ** lCD);

       int updateCellFlux(Tensor<double>& fluxHD,
                          Tensor<double>& fluxCD,
                          const Tensor<weights>& allwgtsHD, double ** lHD,
                          const Tensor<weights>& allwgtsCD, double ** lCD);

       int updateVel_Pause(Tensor<vector<vertex>>& phasevel_vert, 
                           Tensor<vector<vertex>>& phasevel_hori, 
                           Tensor<vector<vertex>>& effvel_vert, 
                           Tensor<vector<vertex>>& effvel_hori, 
                           Tensor<vector<vertex>>& solidvel_vert, 
                           Tensor<vector<vertex>>& solidvel_hori, 
                           const Tensor<weights>& allwgtsHD, double ** lHD,
                           const Tensor<weights>& allwgtsCD, double ** lCD);

       int velocitycamera(Tensor<vertexSet>& phasevel_vert, 
                          Tensor<vertexSet>& phasevel_hori, 
                          Tensor<vertexSet>& effvel_vert, 
                          Tensor<vertexSet>& effvel_hori, 
                          Tensor<vertexSet>& solidvel_vert, 
                          Tensor<vertexSet>& solidvel_hori, 
                          Vec * gCD, Vec * gHD, 
                          int t, double dt, int maxIter, 
                          double tolUzawa);

       int getflux(const Tensor<weights>& allwgtsHD, double ** lHD, 
                   const Tensor<weights>& allwgtsCD, double ** lCD, 
                   double **lfHD, double** lfCD);

       int getflux(const Tensor<vertexSet>& phasevel_vert, 
                   const Tensor<vertexSet>& phasevel_hori, 
                   const Tensor<vertexSet>& effvel_vert, 
                   const Tensor<vertexSet>& effvel_hori, 
                   const Tensor<vertexSet>& solidvel_vert, 
                   const Tensor<vertexSet>& solidvel_hori,
                   const Tensor<weights>& allwgtsHD, double ** lHD,
                   const Tensor<weights>& allwgtsCD, double ** lCD,
                   double **lfHD, double **lfCD);

};

#endif
