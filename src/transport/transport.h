#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "weno_multilevel.h"

class Transport{
    public:
        Transport(const MeshInfo &mi, point_index& cell);
        ~Transport();


    private:

        int blayer_;

        // Starting vertx index for this cell
        point_index cell_;
        int startVertx_[2];

        bool WithinBoundary(int i, int j); 

        index_set Onboundary;
 
        index_set InsideCell; 

        vector<int *> adv_rangex_;
        vector<int *> adv_rangey_;

        vector<int *> diff_rangex_;
        vector<int *> diff_rangey_;

        double AdvectionFlux(const MeshInfo& mi);

        double DiffusionFlux(const MeshInfo& mi);

        double AdvectionFluxOnBoundary(const MeshInfo& mi);

        double DiffusionFluxOnBoundary(const MeshInfo& mi);

};

#endif
