#include "mlwenouse.h"

using namespace MLWENO;

int MLWENOUse::AddMLWENOLevel(const std::string& location,
                              const unordered_set<std::string>& selectLevels,
                              MLWENOPrepare * mlpPtr){

    if (mlpPtr == NULL){
        return -1;
    } else {
        multiLevelReconstruction * mlrins = new multiLevelReconstruction(); 

        mlrins->SelectWenoReconstLevel(selectLevels, (*mlpPtr));

        mlrIns_.push_back(mlrins); 

        AssignMap_.insert(std::make_pair(location, mlrIns_.size()-1));

        return 0;
    }

}

void MLWENOUse::AssignWENOStencils(const int& location,
                                   const std::string& level,
                                   const vector<indice>& newReconstMethod){

    mlrIns_.at(location)->ModifyReconstMethod(level, newReconstMethod);   
}

void MLWENOUse::AssignWENOStencils(const std::string& location,
                                   const std::string& level, 
                                   const vector<indice>& newReconstMethod){
    mlrIns_.at(AssignMap_.at(location))->ModifyReconstMethod(level, newReconstMethod);   
}


double MLWENOUse::Evaluate(const vertex& point,
                           const indice& globalCell,
                           const MeshInfo& mi){

    // AssignInstance will decide which MLWENO reconstruction instance 
    // should be used here.

    //return mlrIns_.at(AssignInstance_(globalCell))->EvaluateMLWENO(mi,point,globalCell); 
    return 0.0;
}

int MLWENOUse::AssignInstance_(const indice& globalCell) const{

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
