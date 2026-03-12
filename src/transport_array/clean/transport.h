#ifndef TRANSPORT_CLEAN_H_
#define TRANSPORT_CLEAN_H_

#include "util.h"
#include "reconstruction.h"
#include "tensorstencilpoly.h"

#include "lagrange_tmp.h"

// Advection related functions
double advfunc(const double& u, const vertex& vel, const vertex& unitnormal);

double ufunc(double u);

double dfdu(const double& u);

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work);

double inflow(const vertex& point, const vector<double>& param);

bool isinflow(const double& flux);

int diriBndry(const vector<vertex>& points, vector<double>& value, int flag);

// Diffusion related functions
double diffusionBndry(const vertex& point, const vector<double>& param);

class TransportVariable{

    public:
        TransportVariable(){};
        ~TransportVariable(){};

        // Global vector holding the solution
        Vec sol; 

        bool diffusion = false;

        // Extra evaluated values
        vector<double> cellgauss;
        vector<double> cellcenter;

        vector<vector<double>> samplingp;

        vector<reconstruction*> my_recon;

        // Create Default (3,2) reconstruction
        int CreateDefaultReconstruction(const MeshInfo& mi);

        // Create a reconstruction
        int CreateReconstruction(const MeshInfo& mi, 
                                 int sizelgx, int sizelgy, int orderlg,
                                 int sizesmx, int sizesmy, int ordersm,
                                 vector<indice>& sten_lg_pre,
                                 vector<indice>& sten_sm_pre);

        // Evaluate reconstruction at a given point
        int Evaluate(const MeshInfo& mi, DM dmu);

        // Get cell flux vector
        int cellflux_all(const MeshInfo& mi, double maxv, DM dmu,
                         const vector<vertex>& edgegaussp,
                         const vector<vertex>& edgevel,
								 bool globalLF, int bndryflag,
                         Vec * influx);

        // Print values at edge gauss points out
        int Print(const MeshInfo& mi, const char * fieldname);

        private:
            vector<tensorstencilpoly> stenlg;
            vector<tensorstencilpoly> stensm;

            vector<double> sigma_lg;
            vector<double> sigma_sm;

            // Auxilliary functions
	         int getNeighbors(int xMaxCell,   int yMaxCell, 
                             int xSize,      int ySize,
                             int i, int j, int& edgepos, int& edgeneg,
                             indice& cellpos, indice& cellneg,
					              bool& onbndry);

            int getUniformEdge(const vertexSet& edge, vertexSet& uniformEdge, const indice& cellid,
									    int xSize, int ySize, int xMaxCell, int yMaxCell);
				
            // Diffusion sampling related variables
			   int degree  = -1;	
            int halfpts = -1;
            int numpts  = -1;

            LagrangeBasisDeriv lagDer;

            // Update reconstruction nonlinear weights
            int UpdateRecon(const MeshInfo& mi, double ** locvals);	

            // Evaluate reconstruction at all the gauss points
            int EvaluateEdge(int i, int j, const MeshInfo& mi, double ** locvals);

            // Evaluate reconstruction at other points
            int EvaluateExtra(int i, int j, const MeshInfo& mi, double ** locvals);

            // Compute advective flux
            double advflux_edge(const vector<vertex>& vel, 
                                const vector<double>& uneg,
                                const vector<double>& upos,
                                const vertexSet& edge,
                                bool localLF, double gLF);
/*
            double advflux_edge(const MeshInfo& mi, 
                                const vector<vertex>& vel, 
                                const vector<double>& u,
                                const vector<double>& f,
                                const vertexSet& edge,
                                bool localLF, double gLF);
*/

            // Compute diffusive flux
            double difflux_edge(const vector<double>& sample,
                                const vertexSet& edge,
										  double len, vertex unitNormal,
										  LagrangeBasisDeriv& lagDer, 
										  double dx, int numPts);

            // =====================================================================

            int ExtractThisEdge(const MeshInfo& mi,
                                indice cellneg, int edgeneg,
										  indice cellpos, int edgepos,
                                int gsize,
                                vector<double>& uneg,
                                vector<double>& upos);

            int EvaluateSamples(const vertexSet& edge,
                                const vertex& unitNormal,
                                double dx,
										  int halfpts,
										  int cellidneg, int cellidpos,
										  double ** locvals,
                                vector<double>& samples);

            int EvaluateSamples(const MeshInfo& mi, 
									     const indice& cellid,
										  const vertexSet& edge,
										  const vertex& unitNormal,
										  int xSize, int ySize, 
										  int xMaxCell, int yMaxCell);

            int EvaluateSamplesEdge(const MeshInfo& mi,
                                    int xSize, int ySize, int offset, 
                                    int xMaxCell, int yMaxCell,
					                     double ** locvals);

            int advflux_edge_all(const MeshInfo& mi, 
                                 int xSize, int ySize, int offset,
                                 int xMaxCell, int yMaxCell,
                                 const vector<vertex>& edgegaussp,
                                 const vector<vertex>& edgevel,
                                 vector<double>& edgeflux,
                                 bool localLF, int bndryflag,
                                 double gLF);

            int difflux_edge_all(const MeshInfo& mi, 
                                 int xSize, int ySize, int offset,
                                 int xMaxCell, int yMaxCell,
                                 const vector<vertex>& edgegaussp,
                                 vector<double>& edgeflux);
};

#endif
