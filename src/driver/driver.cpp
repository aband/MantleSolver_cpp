#include "driver.h"

int Driver::UseWeno(){

    // Allocate memory space for mlweno prepare class

    mlpPtr_ = new MLWENO::MLWENOPrepare();

    return 0;
}

int Driver::AddLevel(const int& m, const int& n){
    if (mlpPtr_ == NULL){
        return -1;
    } else {
        mlpPtr_->AddLevel(mi, m, n);
        return 0;
    }
}

int Driver::PrepareTransport(transportType type){

    switch(type) {
        case adv:
            // Advection only
            break;
        case diff:
            // Diffusion only
            break;
        case adv_diff:
            // Advection-diffusion 
            break;
        case adv_diff_react:
            // Advection-diffusion-reaction
            break;
        default:
            cout << "Not a valid transport ..." << endl;
            break;
    }

    return 0;
}
