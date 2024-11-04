#include "advectiveFlux.h"

inline double LFFlux(const vector<double>& u, 
                     const vector<double>& fu,
                     const double& LF){
    // fu[0] corresponding to fuIn , fu[1] corresponding to fuOut
	 // u[0] corresponding to uIn , u[1] corresponding to uOut

    return 0.5*(fu[0]+fu[1] - LF*(u[1] - u[0]));
}

vector<double> advFlux(const MeshInfo& mi,
                       const MLWENO::MLWENOUse& mluIn,
                       const MLWENO::MLWENOUse& mluOut,
                       const vector<vertex>& edge,
                       const vertex& unitNormal,
                       const double& len,
                       const indice& gCellIn,
                       const indice& gCellOut,
                       const std::string& locIn,
                       const std::string& locOut,
                       const vector<double>& LFparam,
                       const vector<double>& direction,
                       const valarray<double>& gpe){

    assert(LFparam.size() == gpe.size());

    vector<double> work;
    work.resize(gpe.size()); 

    for (int g=0; g<gpe.size(); g++){

        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
        double uIn  = mluIn.Evaluate(mapped, gCellIn, mi, locIn);
        double uOut = mluOut.Evaluate(mapped, gCellOut, mi, locOut);

        // Compute Lax_Friedrich flux
        // Here transport equation is hard coded within
        work.at(g) = LFFlux({uIn, uOut}, 
          {direction.at(g)*uIn, direction.at(g)*uOut}, 
          LFparam.at(g));

    }

    return work;
}

vector<double> advFluxBndry(const MeshInfo& mi,
                            const MLWENO::MLWENOUse& mlu,
                            const vector<vertex>& edge,
                            const vertex& unitNormal,
                            const double& len,
                            const indice& gCell,
                            const std::string& loc,
                            const vector<double>& direction,
                            const vector<double>& LFparam,
                            const valarray<double>& gpe,
                            const int& flag){

    assert(LFparam.size() == gpe.size());

    vector<double> work;
    work.resize(gpe.size()); 

    for (int g=0; g<gpe.size(); g++){

    }

    return work;
}

// =========== Implicit =================================
