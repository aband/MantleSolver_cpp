#include "transport.h"

using namespace Transport {

    void advection::Clear() const {
        boundaryType_.clear();

        for (auto it = advRecon_.begin(); it != advRecon_.end(); it++){
            delete it->second;
        }
    }

    void CreateBoundary_() {


    }


    void diffusion::Clear() const {
        boundaryType_.clear();

        for (auto it = diffRecon_.begin(); it != diffRecon_.end(); it++){
            delete it->second;
        }
    }



}
