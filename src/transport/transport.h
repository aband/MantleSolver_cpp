#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "weno_multilevel.h"

// Define transport per cell
class TransportCell{
    public:
        TransportCell(cosnt MeshInfo& mi, point_index& cell);
        ~TransportCell();

    private:
        // Starting vertx index for this cell
        point_index cellIndex_;
        int startVertx_[2];

        vector<int *>  advRangex_;
        vector<int *>  advRangey_;

        vector<int *>  diffRangex_;
        vector<int *>  diffRangey_;

        WenoReconstruction * advWr_;
        WenoReconstruction * diffWr_;

        // Functions used for flux computation
        double ComputeAdvectionFlux(const MeshInfo& mi);

        double ComputeDiffusionFlux(const MeshInfo& mi);

        double ComputeAdvectionFluxOnBoundary(const MeshInfo& mi);

        double ComputeDiffusionFluxOnBoundary(const MeshInfo& mi);

}

class Transport{
    public:
        Transport(const MeshInfo &mi, point_index& cell);
        ~Transport();

    private:
        int blayer_;

        bool WithinBoundary(int i, int j); 

        index_set Onboundary;

        index_set InsideCell; 

};

#endif
