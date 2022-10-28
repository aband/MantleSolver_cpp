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
      



}

// ========================================================================

Transport::Transport(){


}

Transport::~Transport(){

}


Transport::WithinBoundary(int i, int j){

    if (i<0+blayer_ || i>mi.globalsize[0]-blayer_ || j<0+blayer_ || j>mi.globalsize[1]-blayer_){
        return true;
    } else {
        return false;
    }

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


