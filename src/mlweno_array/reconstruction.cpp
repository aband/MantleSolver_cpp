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
const bool validsten(int ){


}

// =========================================================================

int reconstruction::init(){

    count = sten_lg.size() + sten_sm.size() + use_sten_const; 

    linwgts.clear();
    linwgts.resize(count);

    nonlinwgts.clear();
    nonlinwgts.resize(count);

    return 1;
}

int reconstruction::extractsigma(const vector<double>& sigma_lg, 
                                 const vector<double>& sigma_sm){



    return 1;
}

int reconstruction::setWgts(const vector<double>& stensigma, double h0){

    double sum = 0.0;

    nonlinwgts.clear(); 

    for (int l=0; l<nonlinwgts.size(); l++){
        nonlinwgts.at(l) = linwgts.at(l) / pow(stensigma.at(l) + epsilon*h0*h0, s*stenorder.at(l) + geteta(stenorder.at(l)));
        sum += nonlinwgts.at(l);
    }

    for(int l=0; l<nonlinwgts.size(); l++){
        nonlinwgts.at(l) /= sum;
    }

    return 1;
}
