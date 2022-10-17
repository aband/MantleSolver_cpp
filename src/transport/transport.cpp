/*
 * Define a cpp function for degenerate hyperblic equation
 */

#include "../../include/transport.h"

// Initialize a transport object with respect to a given index number of cell
// Cell index represents global cell index
TransportCell::TransportCell(const MeshInfo& mi, point_index& cellIndex){

    assert(cell.size() == 2);

    cellIndex_ = cell;

    startVertx_[0] = cell[0] + mi.ghost_vertx[0];
    startVertx_[1] = cell[1] + mi.ghost_vertx[1];

}

TransportCell::GetAdvStencil(vector<int *> rangex, vector<int *> rangey){
    advRangex_ = rangex;
    advRangey_ = rangey;


}

TransportCell::GetDiffStencil(vector<int *> rangex, vector<int *> rangey){
    diffRangex_ = rangex;
    diffRangey_ = rangey;


}

TransportCell::PrepareWenoReconstruction(const MeshInfo& mi){

    advWr_ = new WenoReconstruction(mi,advLinWeights,advRangex_,advRangey_,cellIndex_);;
    diffWr_ = new WenoReconstruction(mi,diffLinWeights,diffRangex_,diffRangey_,cellIndex_);;

}

// =========================================================================



Transport::DetermineBoundarylayer(){


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

Transport::SetupAdvFluxReconstructionInfo(const meshinfo& mi, vector<double>& linweights, point_index& target){

    advWr.push_back(new WenoReconstruction(mi,linWeights,advRangex,advRangey,target)); 

}

Transport::SetupDiffFluxReconstructionInfo(const meshinfo& mi, vector<double>& linweights, point_index& target){

    advWr.push_back(new WenoReconstruction(mi,linWeights,diffRangex,diffRangey,target));

}
