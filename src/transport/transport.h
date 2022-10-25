#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "weno_multilevel.h"

enum boundaryType {none = 0, inflow = 1, outflow = 2};

// Define transport per cell
class TransportCell{
    public:
        TransportCell(cosnt MeshInfo& mi, point_index& cell);
        ~TransportCell();

        GetAdvStencil(vector<int>& order);
        GetAdvStencil(vector<int *> rangex, vector<int *> rangey);

        GetDiffStencil(vector<int>& order);
        GetDiffStencil(vector<int *> rangex, vector<int *> rangey);

    private:
        // Starting vertx index for this cell
        point_index cellIndex_;
        int startVertx_[2];

        // Indicate type of boundary
        boundaryType boundaryType_ = none;

        // Define a master stencil
        vector<int *>  advRangex_;
        vector<int *>  advRangey_;

        vector<int *>  diffRangex_;
        vector<int *>  diffRangey_;

        // Define weno stencils for reconstruction
        WenoReconstruction * advWr_;
        WenoReconstruction * diffWr_;

        // Functions used for flux computation
        double ComputeAdvectionFlux(const MeshInfo& mi);

        double ComputeDiffusionFlux(const MeshInfo& mi);
}

class Transport{
    public:
        Transport();
        ~Transport();

    private:
        int blayer_ = 1;

        bool WithinBoundary(int i, int j); 

        index_set Onboundary_;

        index_set InteriorCell_; 

        void FindBoundary(const MeshInfo& mi);

};

#endif
