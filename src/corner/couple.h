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

extern "C"{
#include "mesh.h"
#include "output.h"
//#include "cgns_io.h"
}

// MFEM parameter header file
#include "myFunc.h"

// Coupling require the following values evaluated on the cell edges and cell
// face quadrature points
// This is a new strstegy. Everything evaluated at each gauss points, being 
// stored and passed between flow and transport solver.

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

         Vec globalCD, globalHD;

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

        // Actual coupling functions
        int computePorosity();

        int computePorosity_phase();

        int printedgeporosity(int mark);

        // edge velocity
        int printedgevel(int mark, const vector<vertex>& stokesvel,
                                   const vector<vertex>& darcyvel);

        // Coupling variables
        vector<vertex> edgegauss;
        vector<vertex> cellgauss;
        vector<vertex> cellcenter;

        vector<double> edgeporo;
        vector<double> cellporo;
        vector<double> average_poro;

        vector<vertex> phasevel;
        vector<vertex> effvel;
        vector<vertex> solidvel;
        vector<vertex> liquidvel;

    private:
        // Parameters
        double L_, H_;
        int M_, N_;

        // Primary variables

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
