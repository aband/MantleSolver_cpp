#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "weno_multilevel.h"

enum boundaryType {inflow, outflow};
enum edgeIndex {West, South, East, North};

class TransportCell{
    public:
        TransportCell(const MeshInfo& mi, point_index& cellIndex);
        ~TransportCell();

        // Indicate if the given cell is on boundary
        bool boundaryflag;

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

        // Boundary information stores which edge is on the boundary and 
        // the corresponding boundary type: inflow or outflow
        vector<pair <int, boundaryType>> boundaryInfo_;

        bool withinBoundary_(const MeshInfo& mi); 

        void identifyBoundary_(double * horieffVel, double * vertEffVel);

}

// Define transport per cell
class Transport{
    public:
        Transport(const MeshInfo& mi);
        ~Transport();

    private:
        // vector holding cell index of boundary cells and interior cells
        vector< transportCell *> onboundarycell_;
        vector< transportCell *> interiorcell_;

        // vector holding reconstruction method of advective flux and diffusion flux
        vector<WenoReconstruction *> advWr_;
        vector<WenoReconstruction *> diffwr_;

};

#endif
