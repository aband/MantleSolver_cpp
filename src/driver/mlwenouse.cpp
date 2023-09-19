#include "mlweno.h"

using namespace MLWENO;

void MLWENOUse::AddMLWENOInstance(const unordered_set<std::string>& selectLevels,
                                  MLWENOPrepare * mlpPtr){

    multiLevelReconstruction * mlrins = new multiLevelreconstruction(); 

    mlrins->SelectWenoReconstLevel(selectLevels, (*mlpPtr));

    mlrIns.push_back(mlrins); 

}

int AssignInstance(const indice& globalCell){

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



double Evaluate(const vertex& point,
                const indice& globalCell){

    // AssignInstance will decide which MLWENO reconstruction instance 
    // should be used here.

    return mlrIns_.at(AssignInstance(globalCell))->EvaluateMLWENO(mi,point,globalCell); 
}
