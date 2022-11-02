/*
 * Define a cpp function for degenerate hyperblic equation
 */

#include "../../include/transport.h"

TransportCell::TransportCell(const MeshInfo& mi, point_index& cellIndex){

    assert(cellIndex.size() == mi.ghost_vertx.size());

    // Local cell index without considering ghost layers
    localCellIndex_ = cellIndex;
    localCellIndexFlat_ = localCellIndex_[1]*mi.localsize[0] + localCellIndex_[0];

    // local cell index considering ghost layers
    localCellIndexGhost_[0] = localCellIndex_[0] + mi.ghost_cell[0];
    localCellIndexGhost_[1] = localCellIndex_[1] + mi.ghost_cell[1];
    localCellIndexGhostFlat_ = localCellIndexGhost_[0] * (mi.localsize[0]+2*mi.ghost_cell[0]) +
                               localCellIndexGhost_[1];
    // global cell index
    globalCellIndex_[0] = localCellIndex_[0] + mi.localstart[0]; 
    globalCellIndex_[1] = localCellIndex_[1] + mi.localstart[1]; 
    globalCellIndexFlat_ = globalCellIndex[1]*mi.globalsize[0] + globalCellIndex[0];

    // local edge index
    localNVertiEdge_ = mi.localsize[0] + 1;
    localNHoriEdge_ = mi.localsize[0];

    globalNVertiEdge_ = mi.globalsize[0] + 1;
    globalNHoriEdge_ = mi.globalsize[0];

    // edge indexing order
    // |-3-|
    // 0   2
    // |-1-|

    localEdgeIndex_[1] = localCellIndex_[1]*localNHoriEdge_ + localCellIndex_[0];
    localEdgeIndex_[3] = (localCellIndex_[1]+1)*localNHoriEdge_ + localCellIndex_[0];

    localEdgeIndex_[0] = localCellIndex_[1]*localNVertiEdge_ + localCellIndex_[0] + 1;
    localEdgeIndex_[2] = localCellIndex_[1]*localNVertiEdge_ + localCellIndex_[0];

    globalEdgeIndex_[1] = globalCellIndex_[1]*globalNHoriEdge_ + globalCellIndex_[0];
    globalEdgeIndex_[3] = (globalCellIndex_[1]+1)*globalNHoriEdge_ + globalCellIndex_[0];

    globalEdgeIndex_[0] = globalCellIndex_[1]*globalNVertiEdge_ + globalCellIndex_[0] + 1;
    globalEdgeIndex_[2] = globalCellIndex_[1]*globalNVertiEdge_ + globalCellIndex_[0];

    // Assign boundary information
    if (withinBoundary_()) {

        boundaryflag = true;

        identifyBoundary_(horiEffVel, vertEffVel);

    } else {

        boundaryflag = false;

    }

}

TransportCell::identifyBoundary_(double * horiEffVel, double * vertEffVel){

    // Find which edges are on the boundary
    // Identify boundry types at the same time

    if (globalCellIndex_[0] == 0) {

        edgeIndex ei = West; 
        pair tmp;
        if (horiEffVel[localEdgeIndex_[ei]] > 0){
            tmp = make_pair(ei,outflow);    
        } else {
            tmp = make_pair(ei,inflow);    
        }

        boundaryInfo_.push_back(tmp);

   } else if (globalCellIndex_[0] == mi.globalsize[0]-1) {

        edgeIndex ei = East;
        pair tmp;
        if (horiEffVel[localEdgeIndex_[ei]] > 0){
            tmp = make_pair(ei,outflow);    
        } else {
            tmp = make_pair(ei,inflow);    
        }

        boundaryInfo_.push_back(tmp);

   } else if (globalCellIndex_[1] == 0){

        edgeIndex ei = South;
        pair tmp;
        if (horiEffVel[localEdgeIndex_[ei]] > 0){
            tmp = make_pair(ei,outflow);    
        } else {
            tmp = make_pair(ei,inflow);    
        }

         boundaryInfo_.push_back(tmp);

   } else {

        edgeIndex ei = North;
        pair tmp;
        if (horiEffVel[localEdgeIndex_[ei]] > 0){
            tmp = make_pair(ei,outflow);    
        } else {
            tmp = make_pair(ei,inflow);    
        }
 
        boundaryInfo_.push_back(tmp);

    }

}

TransportCell::withinBoundary_(const MeshInfo& mi){

    if (globalCellIndex_[0]==0 || globalCellIndex_[0]==mi.globalsize[0]-1 || 
        globalCellIndex_[1]==0 || globalCellIndex_[1]==mi.globalsize[1]-1 ){
        return true;
    } else {
        return false;
    }

}

// ========================================================================

Transport::Transport(const MeshInfo& mi){

    for (int j=0; j<mi.localsize[1]; j++){
    for (int i=0; i<mi.localsize[0]; i++){

        point_index cid {i,j};

        transportCell * tmpPtr = new transportCell(mi, cid);

        if (tmpPtr->boundaryflag){
            onboundarycell_.push_back(tmpPtr); 
        } else {
            interiorcell_.push_back(tmpPtr);
        }

    }}

}

Transport::~Transport(){

    // delete vectors of pointers
    for (auto ptr : onboundarycell_){
        delete ptr;
    }
    onboundarycell_.clear();

    for (auto ptr : interiorcell_){
        delete ptr;
    }
    interiorcell_.clear();

    for (auto ptr: advWr_){
        delete ptr; 
    }
    advWr_.clear();

    for (auto ptr: diffWr_){
        delete ptr; 
    }
    diffWr_.clear();

}

Transport::FindBoundary(const MeshInfo& mi){

    for (int j=mi.localsize[1]-mi.ghost_cell[1]; j<mi.localsize[1]+mi.ghost_cell[1];j++){
    for (int i=mi.localsize[0]-mi.ghost_cell[0]; i<mi.localsize[0]+mi.ghost_cell[0];i++){

        valarray<int> currentCell = {i,j};

        if (WithinBoundary(i,j)) {
            Onboundary_.push_back(currentCell);
        } else {
            InteriorCell_.push_back(currentCell); 
        } 

    }}

}

// Separate boundary cells with inner cells
Transport::SeparateBoundary(const MeshInfo& mi){

    for (int j=0; j<mi.localsize[1]; j++){
    for (int i=0; i<mi.localsize[0]; i++){
        int currenti = cell_[0] + i;
        int currentj = cell_[1] + j;
        if (WithinBoundary(currenti, currentj)){
            Onboundary.push_back({i,j});
        } else {
            InteriorCell.push_back({i,j});
        }
 
    }}

}


