#include "reconstruction.h"
#include "mlwenouse.h"

using namespace MLWENO;

void MLWENOUse::AddMLWENOInstance(const unordered_set<std::string>& selectLevels,
                                  MLWENOPrepare * mlpPtr){

    multiLevelReconstruction * mlrins = new multiLevelReconstruction(); 

    mlrins->SelectWenoReconstLevel(selectLevels, (*mlpPtr));

    mlrIns_.push_back(mlrins); 

}

double MLWENOUse::Evaluate(const vertex& point,
                           const indice& globalCell,
                           const MeshInfo& mi){

    // AssignInstance will decide which MLWENO reconstruction instance 
    // should be used here.

    return mlrIns_.at(AssignInstance_(globalCell))->EvaluateMLWENO(mi,point,globalCell); 
}

int MLWENOUse::AssignInstance_(const indice& globalCell){

    // No special treatment on boundary as default 

    return 0;

    // ==== template ========================
    // These boundaryx function can be defined as inline functions in this file.
//    if (boundary1(globalCell)) {
//        return 0;
//    } else if (boundary2(globalCell)){
//        return 1;
//    } else {
//        return 2;
//    }

}
