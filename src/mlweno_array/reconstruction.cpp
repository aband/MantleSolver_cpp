#include "reconstruction.h"

const int geteta(int r){

    if (r==1){
        return 1;
    } else if (r==2){
        return 3;
    } else {
        return 4;
    }
}

// Return whether it is a valid stencil index
// In serial manner
const bool validsten(int sizex, int sizey, 
                     const MeshInfo& mi, indice stencil){

    int M = mi.MPIglobalCellSize[0] - sizex+1;
    int N = mi.MPIglobalCellSize[1] - sizey+1;

    if (stencil[0] < 0 || stencil[1] < 0 || stencil[0] > M-1 || stencil[1] > N-1 ){
        // it is not a valid stencil
        return false;
    } else {
        // it is a valid stencil
        return true;
    }
}

// =========================================================================

int reconstruction::init(int sizex_sm, int sizey_sm,
                         int sizex_lg, int sizey_lg,
                         int order_sm, int order_lg,
								 const vector<indice>& sten_lg_pre,
                         const vector<indice>& sten_sm_pre,
                         const MeshInfo& mi, indice start){

    sten_lg.clear();
    sten_sm.clear();

    for (auto it: sten_lg_pre){
        if (validsten(sizex_lg, sizey_lg, mi, start + it)){
				sten_lg.push_back(it);
				linwgts_lg.push_back(1);
        }
    }

    for (auto it: sten_sm_pre){
        if (validsten(sizex_sm, sizey_sm, mi, start + it)){
				sten_sm.push_back(it);
            linwgts_sm.push_back(1);
        }
    }

    r_lg = order_lg + 1;
    r_sm = order_sm + 1;

    return 1;
}

int reconstruction::extractsigma(const vector<double>& sigma_lg, 
                                 const vector<double>& sigma_sm){



    return 1;
}

int reconstruction::setWgts(const vector<double>& stensigma_lg, 
                            const vector<double>& stensigma_sm,
									 double h0){

    double sum = 0.0;

    nonlinwgts_lg.clear(); 
    nonlinwgts_sm.clear(); 

    for (int l=0; l<nonlinwgts_lg.size(); l++){
        nonlinwgts_lg.at(l) = linwgts_lg.at(l) / pow(stensigma_lg.at(l) + epsilon*h0*h0, s*r_lg + geteta(r_lg));
        sum += nonlinwgts_lg.at(l);
    }

    for (int l=0; l<nonlinwgts_sm.size(); l++){
        nonlinwgts_sm.at(l) = linwgts_sm.at(l) / pow(stensigma_sm.at(l) + epsilon*h0*h0, s*r_sm + geteta(r_sm));
        sum += nonlinwgts_sm.at(l);
    }

    if (use_sten_const){
        nonlinwgts_const = linwgts_const / pow(0.0 + epsilon*h0*h0, s*r_const + geteta(r_const));
    }

    for (int l=0; l<nonlinwgts_lg.size(); l++){
        nonlinwgts_lg.at(l) /= sum;
    }

    for (int l=0; l<nonlinwgts_sm.size(); l++){
        nonlinwgts_sm.at(l) /= sum;
    }

    return 1;
}
