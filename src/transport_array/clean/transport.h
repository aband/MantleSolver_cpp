#ifndef TRANSPORT_CLEAN_H_
#define TRANSPORT_CLEAN_H_

#include "util.h"
#include "reconstruction.h"
#include "tensorstencilpoly.h"

double advfunc(const double& u, const vertex& vel, const vertex& unitnormal);

double ufunc(double u);

double dfdu(const double& u);

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work);

double inflow(const vertex& point, const vector<double>& param);

bool isinflow(const double& flux);

int diriBndry(const vector<vertex>& points, vector<double>& value, int flag);

class TransportVariable{

    public:
        TransportVariable(){};
        ~TransportVariable(){};

        // Global vector holding the solution
        Vec sol; 

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

            vector<reconstruction*> my_recon;

            // Update reconstruction nonlinear weights
            int UpdateRecon(const MeshInfo& mi, double ** locvals);	

            // Evaluate reconstruction at all the gauss points
            int EvaluateEdge(int i, int j, const MeshInfo& mi, double ** locvals);

            // Compute advective flux
            double advflux_edge(const MeshInfo& mi, 
                                const vector<vertex>& vel, 
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


            // =====================================================================

            int ExtractThisEdge(const MeshInfo& mi,
                                indice cellneg, int edgeneg,
										  indice cellpos, int edgepos,
                                int gsize,
                                vector<double>& uneg,
                                vector<double>& upos);

            int advflux_edge_all(const MeshInfo& mi, 
                                 int xSize, int ySize, int offset,
                                 int xMaxCell, int yMaxCell,
                                 const vector<vertex>& edgegaussp,
                                 const vector<vertex>& edgevel,
                                 vector<double>& edgeflux,
                                 bool localLF, int bndryflag,
                                 double gLF);

};

#endif
