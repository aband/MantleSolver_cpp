/*
 * Define a cpp function for degenerate hyperblic equation
 */

#include "../../include/transport.h"

// Initialize a transport object with respect to a given index number of cell
// Cell index represents global cell index
Transport::Transport(point_index& cell){

    assert(cell.size() == 2);

    cell_ = cell;
}

Transport::WithinBoundary(int i, int j){

    if (i<0+blayer_ || i>mi.globalsize[0]-blayer_ || j<0+blayer_ || j>mi.globalsize[1]-blayer_){
        return true;
    } else {
        return false;
    }

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
            InsideCell.push_back({i,j});
        }
 
    }}

}
