#ifndef TRANSPORT_CLEAN_H_
#define TRANSPORT_CLEAN_H_

#include "util.h"
#include "reconstruction.h"
#include "tensorstencilpoly.h"

double advfunc(const double& u, const vertex& vel, const vertex& unitnormal);

double dfdu(const double& u);

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work);

double inflow(const vertex& point, const vector<double>& param);

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

        // Compute advective flux
        double advflux_edge(const MeshInfo& mi, 
                            const vector<vertex>& vel, 
                            const vector<double>& uneg,
                            const vector<double>& upos,
                            const vector<double>& fneg,
                            const vector<double>& fpos,
                            const vertexSet& edge,
                            bool localLF, double gLF);

        double advflux_edge(const MeshInfo& mi, 
                            const vector<vertex>& vel, 
                            const vector<double>& u,
                            const vector<double>& f,
                            const vertexSet& edge,
                            bool localLF, double gLF);

        // Compute diffusive flux

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
};

#endif
