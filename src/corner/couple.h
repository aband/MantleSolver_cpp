#ifndef COUPLE_H_
#define COUPLE_H_

#include <petsc.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include "integral.h"
//#include "eutectic.h"
#include "eutectic_rescaled.h"
#include "input.h"
#include "util.h"

#include "serial_solver.h"
#include "transport.h"

extern "C"{
#include "mesh.h"
#include "output.h"
//#include "cgns_io.h"
}

// MFEM parameter header file
#include "myFunc.h"

// MLWENO part
#include "tensorstencilpoly.h"
#include "reconstruction.h"

// Coupling require the following values evaluated on the cell edges and cell
// face quadrature points
// This is a new strstegy. Everything evaluated at each gauss points, being 
// stored and passed between flow and transport solver.

// Initialization functions for composition and enthalpy
double InitCD(const valarray<double>& point,
              const vector<double>& param);

double InitHD(const valarray<double>& point,
              const vector<double>& param);

// Create output filename 
char * GetFilenameAdd(const char * fieldname, const char * add, int mark);

char * GetFilename(const char * fieldname, int mark);

// Edge values are always vertical edges first then horizontal edges
class couple {

    public:

        couple() {};
        ~couple() {}

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
        Phase myPhase;
        int CreatePhase();
        int ShowPhase();
        int withUnit; 

        // Phase quantities on the edge
        vector<EUTECTIC::PhaseComp> phasequantityedge;

        // cell averaged phase quantities
        vector<EUTECTIC::PhaseComp> phasequantitycell;

        int computePorosityEdgePhase(const MeshInfo& mi, int xSize, int ySize, int offset,
								             int xMaxCell, int yMaxCell, 
												 TransportVariable& H, TransportVariable& C);

        int ExtractThisEdgeGauss(vector<vertex>& thisedgegauss, int gsize, 
								         int i, int j, int xSize, int ySize, int offset);

        int EvaluateThisEdgePhase(int gsize, int i, int j, int xSize, int ySize, int offset, 
								            const vector<double>& Hneg, const vector<double>& Cneg, 
								            const vector<double>& Hpos, const vector<double>& Cpos);

        int computePorosityCellPhase(TransportVariable& H,
												 TransportVariable& C);

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

        int printGaussPoints();

        int printCellGrids();

        /**!
         * Create boundary condition vectors
         */
        int PrepareFlow();

        /**!
         * Create two transport variables 
         */
        int PrepareTransport(TransportVariable& H,
                             TransportVariable& C,
                             double (*initHD)(const valarray<double>& point, 
                                              const vector<double>& param),
                             double (*initCD)(const valarray<double>& point, 
                                              const vector<double>& param));

        /**!
         * Read vector as input
         */
        int ReadVectorTransport(Vec * H, Vec * C, 
                                const char * fileH, 
										  const char * fileC,
										  int mark);

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
         * Split computational domain into different phase regions
         */
        //int PhaseSplit(TransportVariable& H, 
        //               TransportVariable& C);

        // Actual coupling functions
        int computePorosity();

        int computePorosity_phase(TransportVariable& H,
                                  TransportVariable& C);

        int calculatePhaseVel(const vector<vertex>& stokesvel,
								      const vector<vertex>& relativevel);

        int alterporosity(double scale);

        // Adjust edge porosity according to the cell-averaged porosity
        int adjustEdgePorosity();

        int edgegaussToCellIndex(int index, vector<indice>& neighbor);

        int printedgeporosity(int mark);

        int printphase(int mark);

        int printpressure(int mark, Vec * sp, Vec * dp, bool printpressure, 
								const char * spname, const char * dpname);

        // edge velocity
        int printedgevel(int mark, const vector<vertex>& stokesvel,
                                   const vector<vertex>& darcyvel);

        // cell center scalar
        int printCellScalar(Vec * sol, const char * fieldname, int mark);

        int printdivmass(int mark, const char * fieldname, 
								const vector<vertex>& stokesvel,
								const vector<vertex>& darcyvel);

        int printsepmass(int mark, const char * eq1, 
					                    const char * eq2, 
											  const vector<vertex>& stokesvel,
											  const vector<vertex>& darcyvel,
											  const vector<double>& stokesp,
											  const vector<double>& darcyp);

        /**!
         * Coordinates evaluated at different positions
         */
        vector<vertex> edgegauss;
        vector<vertex> cellgauss;
        vector<vertex> cellcenter;

        /**!
         * Porosity evaluated at different positions
         */
        vector<double> edgeporo;
        vector<double> cellporo;
        vector<double> average_poro;

        int examineFullPorosity(int t);

        vector<int> cellphase;
		  vector<double> meltp;
        vector<double> currentTemp;

        int assignTempVec(Vec * temp);

        int setlatentVec(Vec * latent);

        int AssignPorosityVec(Vec * porovec);

        int expandporosity(TransportVariable& poroex, int xMaxCell, int yMaxCell, int xSize, int ySize, int offset);

        int expandporosity();

        /**!
         * Different velocities
         */
        vector<vertex> phasevel;
        vector<vertex> effvel;
        vector<vertex> solidvel;
        vector<vertex> liquidvel;

    private:
        // Parameters
        double L_, H_;
        int M_, N_;

        // scalar values on edges
        int printedgeval(int mark, const vector<double>& val,
                                   const char * fieldname);

        // vector on edges
        int printedgeval(int mark, const vector<vertex>& val,
                                   const char * fieldname);

        // cell center values
        int printcellval(int mark, const vector<double>& val,
                                   const char * fieldname);
};

#endif
