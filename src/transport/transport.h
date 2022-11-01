#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "weno_multilevel.h"

enum boundaryType {none = 0, inflow = 1, outflow = 2};

class TransportCell{
    public:
        TransportCell(const MeshInfo& mi, point_index& cellIndex);
        ~TransportCell();

    private:
        // The prefix local and global are 
        // referring to parallel local and global
        point_index localCellIndex_;
        int localCellIndexFlat_;

        point_index localCellIndexGhost_;
        int localCellIndexFlatGhost_;

        point_index globalCellIndex_;
        int globalCellIndexFlat_;

        // Number of vertical edges and horizontal edges per local row
        int localNVertiEdge_;
        int localNHoriEdge_;

        // Number of vertical edges and horizontal edges per global row
        int globalNVertiEdge_;
        int globalNHoriEdge_;

        int[4] localEdgeIndex_;
        int[4] globalEdgeIndex_;

        // Boundary condition

}

// Define transport per cell
class Transport{
    public:
        Transport();
        ~Transport();

    private:

        vector<WenoReconstruction *> advWr_;
        vector<WenoReconstruction *> diffwr_;

        int blayer_ = 1;

        bool WithinBoundary_(int i, int j); 

        vector< pair<point_index,boundaryType> > ID_;

        void CreateID_(const MeshInfo& mi);
};

#endif
