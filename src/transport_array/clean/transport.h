#ifndef TRANSPORT_CLEAN_H_
#define TRANSPORT_CLEAN_H_

#include "util.h"
#include "reconstruction.h"
#include "tensorstencilpoly.h"

using edgevals = std::vector<std::array<std::array<double, 3>, 4>>;

class TransportVariable{

    public:
        TransportVariable(){};
        ~TransportVariable(){};

        // Reconstructed values on 
        edgevals cv;

        // Global vector holding the solution
        Vec sol; 

        // Create a reconstruction
        int CreateReconstruction(const MeshInfo& mi, 
                                 int sizelgx, int sizelgy, int orderlg,
                                 int sizesmx, int sizesmy, int ordersm,
                                 vector<indice>& sten_lg_pre,
                                 vector<indice>& sten_sm_pre);

        private:
            vector<tensorstencilpoly> stenlg;
            vector<tensorstencilpoly> stensm;

            vector<double> sigma_lg;
            vector<double> sigma_sm;

};

#endif
