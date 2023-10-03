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

// Assign WENO stencils to different mlweno instance pointers
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

// Update Nonlinear weights
void MLWENOUse::UpdateNonLinearWgts(const MeshInfo& mi, 
                                    const int& location,
                                    const std::string& weightType,
                                    bool (*func)(const indice& globalCell,
                                                 const MeshInfo& mi)){
    mlrIns_.at(location)->UpdateNonLinearWgts(mi, weightType, func);
}

void MLWENOUse::UpdateNonLinearWgts(const MeshInfo& mi, 
                                    const std::string& location,
                                    const std::string& weightType,
                                    bool (*func)(const indice& globalCell,
                                                 const MeshInfo& mi)){

    mlrIns_.at(AssignMap_.at(location))->UpdateNonLinearWgts(mi, weightType, func);
}

double MLWENOUse::Evaluate(const vertex& point,
                           const indice& globalCell,
                           const MeshInfo& mi,
                           const int& location){

    return mlrIns_.at(location)->EvaluateMLWENO(mi,point,globalCell); 
}

double MLWENOUse::Evaluate(const vertex& point,
                           const indice& globalCell,
                           const MeshInfo& mi,
                           const std::string& location){

    return mlrIns_.at(AssignMap_.at(location))->EvaluateMLWENO(mi,point,globalCell); 
}

void MLWENOUse::PrintNonLinearWgts(const int& location, 
                                   const MeshInfo& mi){
    mlrIns_.at(location)->PrintNonLinearWgts(mi); 
}

void MLWENOUse::PrintNonLinearWgts(const std::string& location,
                                   const MeshInfo& mi){
    mlrIns_.at(AssignMap_.at(location))->PrintNonLinearWgts(mi);
}

