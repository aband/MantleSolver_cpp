/*
 * Define a cpp function for degenerate hyperblic equation
 */

#include "../../include/transport.h"

Transport::Transport(const MeshInfo& mi){

    for (int j=0; j<mi.localsize[1]; j++){
    for (int i=0; i<mi.localsize[0]; i++){

        point_index cid {i,j};

        transportCell * tmpPtr = new transportCell(mi, cid);

        // Boundary information is stored inside the cell class
        localcells_.push_back(tmpPtr);

    }}

}

Transport::~Transport(){

    // delete vectors of pointers
    // transport cells
    for (auto ptr : localcells_){
        delete ptr;
    }
    localcells_.clear();

    // weno reconstruction class
    for (auto ptr: advWr_){
        delete ptr; 
    }
    advWr_.clear();

    for (auto ptr: diffhoriWr_){
        delete ptr; 
    }
    diffhoriWr_.clear();

    for (auto ptr: diffvertwr_){
        delete ptr;
    }
    diffvertwr_.clear();

    // weno stencils
    for (auto ptr: advRangex_){
        delete ptr;
    }
    advRangex_.clear();

}

void Transport::GetAdvWenoStencils(vector<int *>& rangex, vector<int *>& rangey){

    advRangex_ = rangex;
    advRangex_ = rangey;

}

void Transport::GetDiffHoriWenoStencils(vector<int *>& rangex, vector<int *>& rangey){

    diffHoriRangex_ = rangex;
    diffHoriRangey_ = rangey;
 
}

void Transport::GetDiffVertWenoStencils(vector<int *>& rangex, vector<int *>& rangey){

    diffVertRangex_ = rangex;
    diffVertRangey_ = rangey;

}

vector<double> Transport::SelectAdvLinWeights(const TransportCell& cell){

    vector<double> selectedWeights;
    for (auto s: advStencilSelection_){
        selectedWeights.push_back(s);
    }
    return selectedWeights;
}


void Transport::CreateWenoReconstruction_(const MeshInfo& mi){

    for (auto cell : localcells_){

        vector<double> selectedadvlinweights = SelectAdvLinWeights(cell);
        vector<int *>  selectedrangex     = SelectStencilx(cell);
        vector<int *>  selectedrangey     = SelectStencily(cell);

        WenoReconstruction * tmpPtr = new WenoReconstruction(mi, linWeights, 
                                                             selectedrangex, selectedrangey, 
                                                             cell->GetLocalAdvWenoIndex(mi)); 

        advWr_.push_back(tmpPtr); 
    
    }

}
