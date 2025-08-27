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

extern "C"{
#include "mesh.h"
#include "output.h"
//#include "cgns_io.h"
}

// New transport code
#include "tensorstencilpoly.h"
#include "reconstruction.h"
#include <chrono>
#include "error.h"
#include "advectiveflux.h"
#include "rk.h"

double InitCD(const valarray<double>& point,
              const vector<double>& param);

double InitHD(const valarray<double>& point,
              const vector<double>& param);

double inflow(const valarray<double>& point,
              const vector<double>& param);

int clearall(vector<vertex>& effvel,
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

class Driver {

    public:
        Driver() {};
        ~Driver() {delete myPhase;}

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

         // Two Levels stencils and reconstructions
         int PrepareTransport(double (*initHD)(const valarray<double>& point, 
                                               const vector<double>& param),
                              double (*initCD)(const valarray<double>& point, 
                                               const vector<double>& param));

         vector<tensorstencilpoly> stenlg;
         vector<tensorstencilpoly> stensm;

         // HD and CD will be using the same reconstruction method
         vector<reconstruction> my_recon_HD;  
         vector<reconstruction> my_recon_CD;  

        /**!
         * Create boundary condition vectors
         */
        int PrepareFlow();
 
        /**!
         * Solve flow at the given time step.
         */
        int SolveFlow(int maxIter, double tolUzawa, double ** lHD, double ** lCD);
 
        /**!
         * Scatter distributed vector to all processor.
         * Prepare for velocity reconstruction on gauss points
         */
        int CreateScatterVec();
 
        /**!
         * Assemble reduced linear system.
         * Reducing boundary dofs
         */
        PetscErrorCode ParallelMatrixAssemble(double ** lHD, double ** lCD);

        double HDbottom;
        double CDbottom;

        std::vector<double> parameter;

        int CellAvePorosity(const indice& gcell, double ** lHD, double ** lCD);


        int AssignLocMatStokes(const indice& gcell,
                               double ** lHD,
                               double ** lCD,
                               LocMat * loc);
 
        int AssignLocMatDarcy(const indice& gcell,
                              double ** lHD,
                              double ** lCD,
                              LocMat * loc);
 
        int AssignLocMatCouple(const indice& gcell,
                               double ** lHD,
                               double ** lCD,
                               double& k);

        int computephase(const std::vector<vertex>& gaussp,
                         const vertexSet& edgep,
                         const indice& gcell,
                         double ** lHD, double ** lCD,
                         vector<double>& cs,
                         vector<double>& cl,
                         vector<double>& phi,
                         vector<double>& TD,
                         vector<double>& dTdH,
	                      vector<double>& cd,
                         vector<double>& hd);

        int computeEffVel(const vector<vertex>& gaussp,
                          const vertexSet& edgep,
                          const indice& gcellin, const indice& gcellout,
                          double ** lHD, double ** lCD,
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

        int computeEffVel(const vector<vertex>& gaussp,
                          const vertexSet& edgep,
                          const indice& gcell,
                          double ** lHD,
                          double ** lCD,
                          vector<vertex>& effvel,
                          vector<vertex>& phasevel,
                          vector<vertex>& solidvel,
                          vector<double>& TD,
                          vector<double>& dTdH,
								  vector<double>& CD,
								  vector<double>& HD);


        int computeEffVel_Nonlinear(const vector<vertex>& gaussp,
                                    const vertexSet& edgep,
                                    const indice& gcellin, const indice& gcellout,
                                    double ** lHD, double ** lCD,
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
			          					   vector<double>& HDout,
											   vector<double>& nonlinuin,
											   vector<double>& nonlinuout,
											   vector<double>& nonlinfin,
											   vector<double>& nonlinfout,
											   vector<double>& nonlindfduin,
												vector<double>& nonlindfduout,
												vector<vertex>& nonlinvel);

        // Eat and spit test
        int exactandreconstructTest();

        // Coupled transport functions
        int computeEdgeFlux(vector<double>& edgefluxHD, 
                            vector<double>& edgefluxCD, 
                            double t, 
                            double ** lHD, double ** lCD); 

        int maxIter, tolUzawa, Tmax;
        double dt;

        int getfluxall(Vec * fHD, Vec * fCD, bool updateVel, double t);

        int rk1(); 

        int printGrid();
        int printPhase(bool update, int mark);

        int start = 0;

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

};

#endif
