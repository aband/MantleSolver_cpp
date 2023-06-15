#include "diffusion.h"

// ========== Diffusion ===========================================
// |      Containing functions regarding diffusion flux           |
// ================================================================

/**
 * Return diffusive flux
 */
double SymDiffusion::diffusion::flux_(const double * ru, int n){
    assert(n == 4);
    return ((ru[2]- ru[1])*beta*beta/(2*alpha)-
            (ru[3]- ru[0])*alpha*alpha/(2*beta))/
           (beta*beta-alpha*alpha);
}

double SymDiffusion::diffusion::dflux_(const double * dru, int n){
    assert(n == 4);
    return flux_(dru, n);
}

/**
 * Check if a given target cell is inside the boundary
 * Return an integer that distinguishes different situations near the boundary.
 */
int SymDiffusion::diffusion::Interior_(const int& k, const int& size){

    if (k == 0){
        return 2;
    } else if (k == size){
        return 3;
    } else if (k == 1 || k == size - 1){
        return 1;
    } else {
        return 0;
    }

}

void SymDiffusion::diffusion::UpdateEdgeFlux(const MeshInfo& mi){

    // Update every left and bottom edge for each target cell
    for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] ; j++){
    for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] ; i++){

        indice globalCell {i,j};

        const double scale = sqrt(mi.cellArea.at(FlatIndic(mi,globalCell)));

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, globalCell); 

        // Compute and restore horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};

        edgeHoriFlux_[FlatIndic(mi,globalCell)] = edgeFlux_(mi, globalCell, j, mi.MPIglobalCellSize[1], hori, scale, *(mlrPtrHori_)); 
        
        // Assign additional top edge flux to the physical boundary
        if (j == mi.MPIglobalCellSize[1] - 1){
            hori = {corners.at(2), corners.at(3)};
            indice globalEdge {i,j+1};
            edgeHoriFlux_[FlatIndic(mi,globalEdge)] = -1*edgeFlux_(mi,globalCell, j+1, mi.MPIglobalCellSize[1], hori, scale, *(mlrPtrHori_)); 
        }

        // Compute and restore vertical flux
        vertexSet vert {corners.at(3), corners.at(0)};
        
        edgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalCell)] = edgeFlux_(mi, globalCell, i, mi.MPIglobalCellSize[0], vert, scale, *(mlrPtrVert_));

        // Assign additional right edge flux to the physical boundary
        if (i == mi.MPIglobalCellSize[0] -1 ){
            vert = {corners.at(1), corners.at(2)};
            indice globalEdge {i+1,j};
            edgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalEdge)] = -1*edgeFlux_(mi, globalCell, i+1, mi.MPIglobalCellSize[0], vert, scale, *(mlrPtrVert_));
        }

    }}

}

void SymDiffusion::diffusion::UpdateEdgeFluxDerivative(const MeshInfo& mi){
    // Update every let and bottom edge for each target cell
    for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] ; j++){
    for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] ; i++){
        indice globalCell {i,j};

        const double scale = sqrt(mi.cellArea.at(FlatIndic(mi,globalCell)));

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, globalCell);
    
        // Compute and restore horizontal derivative of flux
        vertexSet hori {corners.at(0), corners.at(1)};

        derivEdgeHoriFlux_[FlatIndic(mi,globalCell)] = derivEdgeFlux_(mi, globalCell, j, mi.MPIglobalCellSize[1], hori, scale, *(mlrPtrHori_));

        if (j == mi.MPIglobalCellSize[1] - 1){
            hori = {corners.at(2), corners.at(3)};
            indice globalEdge {i,j+1};
            unordered_map<int, double> derivedgehori = derivEdgeFlux_(mi, globalCell, j+1, mi.MPIglobalCellSize[1], hori, scale, *(mlrPtrHori_));
            for (auto& d: derivedgehori){
                d.second *= -1;
            }
            derivEdgeHoriFlux_[FlatIndic(mi, globalEdge)] = derivedgehori;
        }

        // Compute and restore vertical flux
        vertexSet vert {corners.at(3), corners.at(0)};
        derivEdgeVertFlux_[FlatIndic(mi,globalCell)] = derivEdgeFlux_(mi, globalCell, j, mi.MPIglobalCellSize[0], vert, scale, *(mlrPtrVert_));

        if (i == mi.MPIglobalCellSize[0] - 1){
            vert = {corners.at(1), corners.at(2)};
            indice globalEdge {i+1,j};
            unordered_map<int, double> derivedgevert = derivEdgeFlux_(mi, globalCell, i+1, mi.MPIglobalCellSize[0], vert, scale, *(mlrPtrVert_));
            for (auto& d:derivedgevert){
                d.second *= -1;
            }
            derivEdgeVertFlux_[FlatIndic(mi,globalEdge)] = derivedgevert;
        }

    }}
}

double SymDiffusion::diffusion::Flux(const MeshInfo& mi, const indice& global){

    double work = 0.0;

    double area = mi.cellArea.at(FlatIndic(mi,global));

    work += edgeHoriFlux_[FlatIndic(mi,global)]; 

    work += -1 * edgeHoriFlux_[FlatIndic(mi,global+mi.faceNormal[2])];

    work += -1 * edgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global+mi.faceNormal[1])];

    work += edgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global)];

    work /= area;
  
    return work;
}

unordered_map<int, double> SymDiffusion::diffusion::derivFlux(const MeshInfo& mi, const indice& global){

    unordered_map<int,double> work;

    // bottom horizontal edge
    unordered_map<int,double> tmp = derivEdgeHoriFlux_[FlatIndic(mi,global)];

    double area = mi.cellArea.at(FlatIndic(mi,global));

    for (auto & derivf : tmp){
        if (work.count(derivf.first) > 0){
            work[derivf.first] += derivf.second/area;
        } else {
            work.insert(std::pair<int,double> (derivf.first, derivf.second/area));
        }
    }

    // top horizontal edge
    tmp = derivEdgeHoriFlux_[FlatIndic(mi,global+mi.faceNormal[2])];

    for (auto & derivf : tmp){
        if (work.count(derivf.first) > 0){
            work[derivf.first] += -1*derivf.second/area;
        } else {
            work.insert(std::pair<int,double> (derivf.first, -1*derivf.second/area));
        }
    }

    // right vertical edge
    tmp = derivEdgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global+mi.faceNormal[1])];

    for (auto & derivf : tmp){
        if (work.count(derivf.first) > 0){
            work[derivf.first] += -1*derivf.second/area;
        } else {
            work.insert(std::pair<int,double> (derivf.first, -1*derivf.second/area));
        }
    }

    // left vertical edge
    tmp = derivEdgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global)];
    for (auto & derivf : tmp){
        if (work.count(derivf.first) > 0){
            work[derivf.first] += derivf.second/area;
        } else {
            work.insert(std::pair<int,double> (derivf.first, derivf.second/area));
        }
    }

    return work;
}

double SymDiffusion::diffusion::boundaryCondition_(double * ru, int n, const int& flag){

    // Reflecive boundary condition
    if (flag == 2) {
        ru[0] = -ru[3];
        ru[1] = -ru[2];
        return flux_(ru,n); 
    } else {
        ru[2] = -ru[1];
        ru[3] = -ru[0];
        return flux_(ru,n);
    }

    // Fixing flux here
    //return fluxvalue...!!!
}

// Compute edge flux using symmetrical one stencil
// Follow numerical scheme mentioned in the old paper
double SymDiffusion::diffusion::edgeFlux_(const MeshInfo& mi,
                            const indice& globalCell, 
                            const int& k,
                            const int& size,
                            const vertexSet& edge,
                            const double& scale,
                            const MLWENO::multiLevelReconstruction& mlrPtr){
    double work = 0.0;

    // Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Find Interpolation positions
    double len = length(edge);
    // Compute unit normal vector pointing outside.
    vertex unitNormal = UnitNormal(edge,len);

    // Two different situations close to the boundary
    switch (Interior_(k,size)){
        case 0 :
            for (int g = 0; g<gpe.size(); g++){
                vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

                double ru[4];
                vertex test = mapped + unitNormal*beta*scale;
                ru[0] = mlrPtr.EvaluateMLWENO(mi,mapped+unitNormal*beta*scale, globalCell);
                ru[1] = mlrPtr.EvaluateMLWENO(mi,mapped+unitNormal*alpha*scale,globalCell);
                ru[2] = mlrPtr.EvaluateMLWENO(mi,mapped-unitNormal*alpha*scale,globalCell);
                ru[3] = mlrPtr.EvaluateMLWENO(mi,mapped-unitNormal*beta*scale, globalCell);

                work += gwe[g] * flux_(ru, 4) * len/2.0;
    
            }

            break;

        case 1 :
            for (int g = 0; g<gpe.size(); g++){
                vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

                double ru[4];

                ru[0] = mlrPtr.EvaluateMLWENO(mi,mapped+unitNormal*scale, globalCell);
                ru[1] = mlrPtr.EvaluateMLWENO(mi,mapped+unitNormal*alpha/beta*scale,globalCell);
                ru[2] = mlrPtr.EvaluateMLWENO(mi,mapped-unitNormal*alpha/beta*scale,globalCell);
                ru[3] = mlrPtr.EvaluateMLWENO(mi,mapped-unitNormal*scale, globalCell);

                work += gwe[g] * flux_(ru, 4) * len/2.0;
    
            }

            break;

        case 2:
            for (int g = 0; g<gpe.size(); g++){
                vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

                double ru[4];

                ru[2] = mlrPtr.EvaluateMLWENO(mi,mapped-unitNormal*alpha*scale,globalCell);
                ru[3] = mlrPtr.EvaluateMLWENO(mi,mapped-unitNormal*beta*scale, globalCell);

                work += gwe[g] * boundaryCondition_(ru, 4, 2) * len/2.0;
    
            }

            break;

        case 3:
            for (int g = 0; g<gpe.size(); g++){
                vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

                double ru[4];

                ru[0] = -1*mlrPtr.EvaluateMLWENO(mi,mapped-unitNormal*beta*scale, globalCell);
                ru[1] = -1*mlrPtr.EvaluateMLWENO(mi,mapped-unitNormal*alpha*scale,globalCell);

                work += gwe[g] * boundaryCondition_(ru, 4, 3) * len/2.0;
    
            }

            break;

        default : 
            cout << "Boundary value not assigned correctly ... " << endl;
            break;

    }

    return work;
}

unordered_map<int, double> SymDiffusion::diffusion::derivEdgeFlux_(const MeshInfo& mi,
                                                     const indice& globalCell,
                                                     const int& k,
                                                     const int& size,
                                                     const vertexSet& edge,
                                                     const double& scale,
                                                     const MLWENO::multiLevelReconstruction& mlrPtr){

    unordered_map<int, double> work;

    // Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Find Interpolation positions
    double len = length(edge);
    // Compute unit normal vector pointing outside.
    vertex unitNormal = UnitNormal(edge,len);

    switch (Interior_(k,size)){
        case 0 :
            for (int g = 0; g<gpe.size(); g++){
                vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
                // all maps should have same set of keys
            }
            break;

        default : 
            cout << "Boundary value not assigned correctly ... " << endl;
            break;
    }

    return work;
}


void SymDiffusion::diffusion::CreateMLWENO(const MLWENO::MLWENOPrepare& mlp, const MeshInfo& mi){
    assert(reconstMethodsVert.empty() == 0);
    assert(reconstMethodsHori.empty() == 0);

    // Create weno levels from reconst method
    // Horizontal edges
    for (const auto& rm: reconstMethodsHori){
        wenoLevelsHori.insert(rm.first);
    }

    mlrPtrHori_->SelectWenoReconstLevel(wenoLevelsHori,mlp); 

    for (const auto& rm: reconstMethodsHori){
        mlrPtrHori_->ModifyReconstMethod(rm.first,rm.second);
    }

    //mlrPtrHori_->SeparateBoundaryLayer(mi,2,{"(1,1)","(2,2)"});
    mlrPtrHori_->SeparateBoundaryLayer(mi);

    // Vertical edges
    for (const auto& rm: reconstMethodsVert){
        wenoLevelsVert.insert(rm.first);
    }

    mlrPtrVert_->SelectWenoReconstLevel(wenoLevelsVert,mlp); 

    for (const auto& rm: reconstMethodsVert){
        mlrPtrVert_->ModifyReconstMethod(rm.first,rm.second);
    }

    //mlrPtrVert_->SeparateBoundaryLayer(mi,2,{"(1,1)","(2,2)"});
    mlrPtrVert_->SeparateBoundaryLayer(mi);
}

void SymDiffusion::diffusion::UpdateNonLinearWgts(const MeshInfo& mi){
    mlrPtrHori_->UpdateNonLinearWgts(mi,2);
    mlrPtrVert_->UpdateNonLinearWgts(mi,2);
}

void SymDiffusion::diffusion::GetInfo(const MeshInfo& mi){

    cout << "MLWENO reconstruction for Horizontal edges ..." << endl;
    mlrPtrHori_->GetInfo();
    mlrPtrHori_->PrintSmoothnessIndicator(mi);
    mlrPtrHori_->PrintNonLinearWgts(mi);

    cout << "MLWENO reconstruction for Vertical edges ..." << endl;
    mlrPtrVert_->GetInfo();
    mlrPtrVert_->PrintSmoothnessIndicator(mi);
    mlrPtrVert_->PrintNonLinearWgts(mi);
}

// =============================================================
// New non symmetric weno reconstruction for diffusion flux
// The new method share the same function name with the 
// previous one, but in different namespace.
// =============================================================

// The computation of flux across the edge is the same as
// symmetrical reconstruction scheme.
double NonSymDiffusion::diffusion::flux_(const std::array<double,4>& ru, 
                                         const double& scale){

    return ( (diffFunc(ru[2])- diffFunc(ru[1]))*beta*beta/(2*alpha)-
             (diffFunc(ru[3])- diffFunc(ru[0]))*alpha*alpha/(2*beta) )/
           (beta*beta-alpha*alpha) / scale;
}

unordered_map<int,double> NonSymDiffusion::diffusion::dflux_(const std::array<unordered_map<int,double>,4>& dru,
                                                             const double& scale){

    derivative work;

    std::array<double,4> coeff;

    for (int k=0; k<4; k++){
        coeff[k] = derivCoeff_[k]/scale;
    }

    for (int k=0; k<4; k++){
        unordered_map_arithmetic(work,dru[k],std::plus<double>(),
                                    coeff[k],std::multiplies<double>());
    }

    return work;
}

unordered_map<int,double> NonSymDiffusion::diffusion::dflux_(const std::array<unordered_map<int,double>,4>& dru,
                                                             const std::array<double,4>& ru,
                                                             const double& scale){

    derivative work;

    std::array<double,4> coeff;

    for (int k = 0; k<4; k++){
        coeff[k] = derivCoeff_[k]*dDiffFunc(ru[k])/scale;
    }

    for (int k=0; k<4; k++){
        unordered_map_arithmetic(work,dru[k],std::plus<double>(),
                                    coeff[k],std::multiplies<double>());
    }

    return work;
}

double NonSymDiffusion::diffusion::boundaryCondition_(std::array<double,4>& ru,
                                                      const std::array<int,2>& posOut,
                                                      const double& scale){

    // Reflective boundary condition
    ru[posOut[0]] = -1*ru[3-posOut[0]];
    ru[posOut[1]] = -1*ru[3-posOut[1]];

    return flux_(ru,scale); 
}

derivative NonSymDiffusion::diffusion::boundaryCondition_(std::array<unordered_map<int,double>,4>& dru,
                                                          std::array<double,4>& ru,
                                                          const std::array<int,2>& posOut,
                                                          const double& scale){

    // Reflective boundary condition
    ru[posOut[0]] = -1*ru[3-posOut[0]];
    ru[posOut[1]] = -1*ru[3-posOut[1]];

    for (int k=0; k<2; k++){
        for (auto & it : dru[3-posOut[0]]){
            dru[posOut[0]].insert(std::pair<int,double>(it.first,it.second*-1));
        }
    }

    return dflux_(dru,ru,scale);
}

/**
 * Check if a given target cell is inside the boundary
 * It is not the same as the old scheme
 * No need to distinguish secondary boundary edge.
 * In other word, it is more resemble a advection interior function.
 */
bool NonSymDiffusion::diffusion::Interior_(const MeshInfo& mi, const indice& target){

    if (target[0] < 0 || target[0] > mi.MPIglobalCellSize[0] - 1 || 
        target[1] < 0 || target[1] > mi.MPIglobalCellSize[1] - 1){
        return false;
    } else {
        return true;
    }
}

// Compute edge flux using non symmetrical stencil
// In this scheme, we are using two different multilevel reonstruction
// for in and out cells.
double NonSymDiffusion::diffusion::edgeFlux_(const MeshInfo& mi,
                                             const indice& globalIn,
                                             const indice& globalOut,
                                             const vertexSet& edge,
                                             const double& scale,
                                             const std::array<int,2>& posIn,
                                             const std::array<int,2>& posOut,
                                             const MLWENO::multiLevelReconstruction& mlrPtrIn,
                                             const MLWENO::multiLevelReconstruction& mlrPtrOut){
    double work = 0.0;

    // Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Find Interpolation positions
    double len = length(edge);
    // Compute unit normal vector pointing outside.
    vertex unitNormal = UnitNormal(edge,len);

    // With new reconstruction scheme
    // We only need to distinguish whether outside cell 
    // is out of the boundary or not.
    if (Interior_(mi,globalOut)){
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
            std::array<double,4> ru;
            ru[posIn[0]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*alpha*scale, globalIn);
            ru[posIn[1]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*beta *scale, globalIn);
            ru[posOut[0]] = mlrPtrOut.EvaluateMLWENO(mi, mapped+unitNormal*beta *scale, globalOut);
            ru[posOut[1]] = mlrPtrOut.EvaluateMLWENO(mi, mapped+unitNormal*alpha*scale, globalOut);
            work += gwe[g] * flux_(ru,scale) * len/2.0;
        } 
    } else {
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
            std::array<double,4> ru;
            ru[posIn[0]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*alpha*scale, globalIn);
            ru[posIn[1]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*beta *scale, globalIn);
            work += gwe[g] * boundaryCondition_(ru,posOut,scale) * len/2.0;
        }
    }

    return work;
}

unordered_map<int, double> NonSymDiffusion::diffusion::derivEdgeFlux_(const MeshInfo& mi,
                                                                      const indice& globalIn,
                                                                      const indice& globalOut,
                                                                      const vertexSet& edge,
                                                                      const double& scale,
                                                                      const std::array<int,2>& posIn,
                                                                      const std::array<int,2>& posOut,
                                                                      const MLWENO::multiLevelReconstruction& mlrPtrIn,
                                                                      const MLWENO::multiLevelReconstruction& mlrPtrOut){
    derivative work;

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);

    vertex unitNormal = UnitNormal(edge, len);

    if (Interior_(mi,globalOut)){
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

            std::array<unordered_map<int,double>, 4> dru;

            std::array<double,4> ru;

            dru[posIn[0]] = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*alpha*scale,globalIn);
            dru[posIn[1]] = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*beta *scale,globalIn);
            dru[posOut[0]] = mlrPtrOut.EvaluateDerivMLWENO(mi,mapped+unitNormal*beta *scale,globalOut);
            dru[posOut[1]] = mlrPtrOut.EvaluateDerivMLWENO(mi,mapped+unitNormal*alpha*scale,globalOut);

            ru[posIn[0]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*alpha*scale, globalIn);
            ru[posIn[1]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*beta *scale, globalIn);
            ru[posOut[0]] = mlrPtrOut.EvaluateMLWENO(mi, mapped+unitNormal*beta *scale, globalOut);
            ru[posOut[1]] = mlrPtrOut.EvaluateMLWENO(mi, mapped+unitNormal*alpha*scale, globalOut);
 
            derivative derivflux = dflux_(dru, ru, scale);

            double multi = gwe[g]*len/2.0;

            unordered_map_arithmetic(work, derivflux, std::plus<double>(), 
                                           multi, std::multiplies<double>());
        }
    } else {
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

            std::array<unordered_map<int,double>, 4> dru;

            std::array<double,4> ru;

            dru[posIn[0]] = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*alpha*scale,globalIn);
            dru[posIn[1]] = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*beta *scale,globalIn);

            ru[posIn[0]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*alpha*scale, globalIn);
            ru[posIn[1]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*beta *scale, globalIn);

            derivative derivflux = boundaryCondition_(dru,ru,posOut,scale);

            double multi = gwe[g]*len/2.0;

            unordered_map_arithmetic(work, derivflux, std::plus<double>(), 
                                           multi, std::multiplies<double>());
 
        }
    }

    return work;
}

// Fully differentiate everything
unordered_map<int, double> NonSymDiffusion::diffusion::derivEdgeFluxFull_(const MeshInfo& mi,
                                                                          const indice& globalIn,
                                                                          const indice& globalOut,
                                                                          const vertexSet& edge,
                                                                          const double& scale,
                                                                          const std::array<int,2>& posIn,
                                                                          const std::array<int,2>& posOut,
                                                                          const MLWENO::multiLevelReconstruction& mlrPtrIn,
                                                                          const MLWENO::multiLevelReconstruction& mlrPtrOut){
    unordered_map<int, double> work;

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);

    vertex unitNormal = UnitNormal(edge, len);

    derivative tmp1, tmp2;

    if (Interior_(mi,globalOut)){
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

            std::array<unordered_map<int,double>, 4> dru;

            std::array<double,4> ru;

            //dru[posIn[0]] = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*alpha*scale,globalIn);
            //dru[posIn[1]] = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*beta *scale,globalIn);
            //dru[posOut[0]] = mlrPtrOut.EvaluateDerivMLWENO(mi,mapped+unitNormal*beta *scale,globalOut);
            //dru[posOut[1]] = mlrPtrOut.EvaluateDerivMLWENO(mi,mapped+unitNormal*alpha*scale,globalOut);

            tmp1 = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*alpha*scale,globalIn);
            tmp2 = mlrPtrIn.EvaluateDerivMLWENOAdd(mi,mapped-unitNormal*alpha*scale,globalIn);

            unordered_map_arithmetic(tmp1, tmp2, std::plus<double>());
            unordered_map_arithmetic(dru[posIn[0]],tmp1,std::plus<double>());

            tmp1.clear(); tmp2.clear();

            tmp1 = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*beta *scale,globalIn);
            tmp2 = mlrPtrIn.EvaluateDerivMLWENOAdd(mi,mapped-unitNormal*beta *scale,globalIn);

            unordered_map_arithmetic(tmp1, tmp2, std::plus<double>());
            unordered_map_arithmetic(dru[posIn[1]],tmp1,std::plus<double>());

            tmp1.clear(); tmp2.clear();


            tmp1 = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped+unitNormal*beta *scale,globalIn);
            tmp2 = mlrPtrIn.EvaluateDerivMLWENOAdd(mi,mapped+unitNormal*beta *scale,globalIn);

            unordered_map_arithmetic(tmp1, tmp2, std::plus<double>());
            unordered_map_arithmetic(dru[posOut[0]],tmp1,std::plus<double>());

            tmp1.clear(); tmp2.clear();

            tmp1 = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped+unitNormal*alpha*scale,globalIn);
            tmp2 = mlrPtrIn.EvaluateDerivMLWENOAdd(mi,mapped+unitNormal*alpha*scale,globalIn);

            unordered_map_arithmetic(tmp1, tmp2, std::plus<double>());
            unordered_map_arithmetic(dru[posOut[1]],tmp1,std::plus<double>());

            tmp1.clear(); tmp2.clear();

            // ==========================================================================

            ru[posIn[0]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*alpha*scale, globalIn);
            ru[posIn[1]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*beta *scale, globalIn);
            ru[posOut[0]] = mlrPtrOut.EvaluateMLWENO(mi, mapped+unitNormal*beta *scale, globalOut);
            ru[posOut[1]] = mlrPtrOut.EvaluateMLWENO(mi, mapped+unitNormal*alpha*scale, globalOut);
 
            unordered_map<int,double> derivflux = dflux_(dru, ru, scale);

            for (auto & derivf : derivflux){
                if (work.count(derivf.first) > 0) {
                    work[derivf.first] += gwe[g]*derivf.second*len/2.0;
                } else {
                    work.insert(std::pair<int, double> (derivf.first, gwe[g]*derivf.second*len/2.0));
                }
            }
        }
    } else {
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

            std::array<unordered_map<int,double>, 4> dru;

            std::array<double,4> ru;

            //dru[posIn[0]] = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*alpha*scale,globalIn);
            //dru[posIn[1]] = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*beta *scale,globalIn);

            tmp1 = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*alpha*scale,globalIn);
            tmp2 = mlrPtrIn.EvaluateDerivMLWENOAdd(mi,mapped-unitNormal*alpha*scale,globalIn);

            unordered_map_arithmetic(tmp1, tmp2, std::plus<double>());
            unordered_map_arithmetic(dru[posIn[0]],tmp1,std::plus<double>());

            tmp1.clear(); tmp2.clear();


            tmp1 = mlrPtrIn.EvaluateDerivMLWENO(mi,mapped-unitNormal*beta *scale,globalIn);
            tmp2 = mlrPtrIn.EvaluateDerivMLWENOAdd(mi,mapped-unitNormal*beta *scale,globalIn);

            unordered_map_arithmetic(tmp1, tmp2, std::plus<double>());
            unordered_map_arithmetic(dru[posIn[1]],tmp1,std::plus<double>());

            // ========================

            ru[posIn[0]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*alpha*scale, globalIn);
            ru[posIn[1]]  = mlrPtrIn.EvaluateMLWENO(mi, mapped-unitNormal*beta *scale, globalIn);

            unordered_map<int,double> derivflux = boundaryCondition_(dru,ru,posOut,scale);

            for (auto & derivf : derivflux){
                if (work.count(derivf.first) > 0) {
                    work[derivf.first] += gwe[g]*derivf.second*len/2.0;
                } else {
                    work.insert(std::pair<int, double> (derivf.first, gwe[g]*derivf.second*len/2.0));
                }
            }
        }
    }

    return work;
}

void NonSymDiffusion::diffusion::UpdateEdgeFlux(const MeshInfo& mi){

    std::array<int,2> posIn;
    std::array<int,2> posOut;

    indice globalOut;

    // Update every left and bottom edge for each target cell
    for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] ; j++){
    for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] ; i++){

        indice globalIn {i,j};

        const double scale = sqrt(mi.cellArea.at(FlatIndic(mi,globalIn)));

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, globalIn); 

        // Compute and restore horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};
        vertexSet vert {corners.at(3), corners.at(0)};

        posIn[0] = 2;posOut[0] = 0;
        posIn[1] = 3;posOut[1] = 1;

        globalOut = {i,j-1};
        // Bottom hoizontal flux
        edgeHoriFlux_[FlatIndic(mi,globalIn)] = edgeFlux_(mi, globalIn, globalOut, hori, scale, posIn, posOut,
                                                          *(mlrPtrHoriUp_), *(mlrPtrHoriDown_));
        globalOut = {i-1,j};
        // Left vertical flux
        edgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalIn)] = edgeFlux_(mi, globalIn, globalOut, vert, scale, posIn, posOut, *(mlrPtrVertRight_), *(mlrPtrVertLeft_));

        if (j==mi.MPIglobalCellSize[1]-1){
            hori = {corners.at(2), corners.at(3)};
            globalOut = {i,j+1};
            posIn[0] = 0; posOut[0] = 2;
            posIn[1] = 1; posOut[1] = 3;
            edgeHoriFlux_[FlatIndic(mi,globalOut)] = edgeFlux_(mi, globalIn, globalOut, hori, scale, posIn, posOut,
                                                                *(mlrPtrHoriDown_), *(mlrPtrHoriUp_));       
        }

        if (i==mi.MPIglobalCellSize[0]-1){
            vert = {corners.at(1), corners.at(2)};
            globalOut = {i+1,j};
            posIn[0] = 0; posOut[0] = 2;
            posIn[1] = 1; posOut[1] = 3;
            edgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalOut)] = edgeFlux_(mi, globalIn, globalOut, vert, scale, posIn, posOut, *(mlrPtrVertLeft_), *(mlrPtrVertRight_));      
        }
    }}
}

void NonSymDiffusion::diffusion::UpdateEdgeFluxDerivative(const MeshInfo& mi){

    std::array<int,2> posIn;
    std::array<int,2> posOut;

    indice globalOut;

    // Update every left and bottom edge for each target cell
    for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] ; j++){
    for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] ; i++){
        indice globalIn {i,j};

        const double scale = sqrt(mi.cellArea.at(FlatIndic(mi,globalIn)));

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, globalIn); 

        // Compute and restore horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};
        vertexSet vert {corners.at(3), corners.at(0)};

        posIn[0] = 2;posOut[0] = 0;
        posIn[1] = 3;posOut[1] = 1;


        globalOut = {i,j-1};
        // Bottom hoizontal flux
        derivEdgeHoriFlux_[FlatIndic(mi, globalIn)] = derivEdgeFlux_(mi,globalIn,globalOut,hori,scale,posIn,posOut,
                                                                     *(mlrPtrHoriUp_), *(mlrPtrHoriDown_));

        globalOut = {i-1,j};
        // Left vertical flux
        derivEdgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalIn)] = derivEdgeFlux_(mi, globalIn, globalOut, vert, scale, posIn, posOut, *(mlrPtrVertRight_), *(mlrPtrVertLeft_));

        if (j==mi.MPIglobalCellSize[1]-1){
            hori = {corners.at(2), corners.at(3)};
            globalOut = {i,j+1};
            posIn[0] = 0; posOut[0] = 2;
            posIn[1] = 1; posOut[1] = 3;
            derivEdgeHoriFlux_[FlatIndic(mi,globalOut)] = derivEdgeFlux_(mi, globalIn, globalOut, hori, scale, posIn, posOut, *(mlrPtrHoriDown_), *(mlrPtrHoriUp_));       
        }

        if (i==mi.MPIglobalCellSize[0]-1){
            vert = {corners.at(1), corners.at(2)};
            globalOut = {i+1,j};
            posIn[0] = 0; posOut[0] = 2;
            posIn[1] = 1; posOut[1] = 3;
            derivEdgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalOut)] = derivEdgeFlux_(mi, globalIn, globalOut, vert, scale, posIn, posOut, *(mlrPtrVertLeft_), *(mlrPtrVertRight_));      
        }
    }}
}

// Full derivative including differentiate non linear weights
void NonSymDiffusion::diffusion::UpdateEdgeFluxDerivativeFull(const MeshInfo& mi){

    std::array<int,2> posIn;
    std::array<int,2> posOut;

    indice globalOut;

    // Update every left and bottom edge for each target cell
    for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] ; j++){
    for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] ; i++){
        indice globalIn {i,j};

        const double scale = sqrt(mi.cellArea.at(FlatIndic(mi,globalIn)));

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, globalIn); 

        // Compute and restore horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};
        vertexSet vert {corners.at(3), corners.at(0)};

        posIn[0] = 2;posOut[0] = 0;
        posIn[1] = 3;posOut[1] = 1;


        globalOut = {i,j-1};
        // Bottom hoizontal flux
        derivEdgeHoriFlux_[FlatIndic(mi, globalIn)] = derivEdgeFluxFull_(mi,globalIn,globalOut,hori,scale,posIn,posOut,
                                                                     *(mlrPtrHoriUp_), *(mlrPtrHoriDown_));

        globalOut = {i-1,j};
        // Left vertical flux
        derivEdgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalIn)] = derivEdgeFluxFull_(mi, globalIn, globalOut, vert, scale, posIn, posOut, *(mlrPtrVertRight_), *(mlrPtrVertLeft_));

        if (j==mi.MPIglobalCellSize[1]-1){
            hori = {corners.at(2), corners.at(3)};
            globalOut = {i,j+1};
            posIn[0] = 0; posOut[0] = 2;
            posIn[1] = 1; posOut[1] = 3;
            derivEdgeHoriFlux_[FlatIndic(mi,globalOut)] = derivEdgeFluxFull_(mi, globalIn, globalOut, hori, scale, posIn, posOut, *(mlrPtrHoriDown_), *(mlrPtrHoriUp_));       
        }

        if (i==mi.MPIglobalCellSize[0]-1){
            vert = {corners.at(1), corners.at(2)};
            globalOut = {i+1,j};
            posIn[0] = 0; posOut[0] = 2;
            posIn[1] = 1; posOut[1] = 3;
            derivEdgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalOut)] = derivEdgeFluxFull_(mi, globalIn, globalOut, vert, scale, posIn, posOut, *(mlrPtrVertLeft_), *(mlrPtrVertRight_));      
        }
    }}
}

double NonSymDiffusion::diffusion::Flux(const MeshInfo& mi, const indice& global){

    double work = 0.0;

    double area = mi.cellArea.at(FlatIndic(mi,global));

    work += edgeHoriFlux_[FlatIndic(mi,global)]; 

    work += -1 * edgeHoriFlux_[FlatIndic(mi,global+mi.faceNormal[2])];

    work += -1 * edgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global+mi.faceNormal[1])];

    work += edgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global)];

    work /= area;

    work *= -1;

    return work;
}

unordered_map<int, double> NonSymDiffusion::diffusion::derivFlux(const MeshInfo& mi, const indice& global){

    derivative work;

    // bottom horizontal edge
    derivative tmp = derivEdgeHoriFlux_[FlatIndic(mi,global)];

    double area = mi.cellArea.at(FlatIndic(mi,global));

    double multi = -1.0/area;

    unordered_map_arithmetic(work,tmp,std::plus<double>(),
                                  multi,std::multiplies<double>()); 

    // top horizontal edge
    tmp = derivEdgeHoriFlux_[FlatIndic(mi,global+mi.faceNormal[2])];

    unordered_map_arithmetic(work,tmp,std::plus<double>(),
                                  area,std::divides<double>());

    // right vertical edge
    tmp = derivEdgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global+mi.faceNormal[1])];

    unordered_map_arithmetic(work,tmp,std::plus<double>(),
                                  area,std::divides<double>());

    // left vertical edge
    tmp = derivEdgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global)];

    unordered_map_arithmetic(work,tmp,std::plus<double>(),
                                  multi,std::multiplies<double>()); 

    return work;
}

void NonSymDiffusion::diffusion::CreateMLWENO(const MLWENO::MLWENOPrepare& mlp, const MeshInfo& mi){
    assert(reconstMethodsVertRight.empty() == 0);
    assert(reconstMethodsVertLeft.empty()  == 0);
    assert(reconstMethodsHoriUp.empty()    == 0);
    assert(reconstMethodsHoriDown.empty()  == 0);

    // Create weno levels form reconst method
    // Vertical right method ===========================================
    for (const auto& rm: reconstMethodsVertRight){
        wenoLevelsVertRight.insert(rm.first);
    }

    mlrPtrVertRight_->SelectWenoReconstLevel(wenoLevelsVertRight, mlp);

    for (const auto& rm: reconstMethodsVertRight){
        mlrPtrVertRight_->ModifyReconstMethod(rm.first,rm.second);
    }

    mlrPtrVertRight_->SeparateBoundaryLayer(mi);

    // Vertical left method ============================================
    for (const auto& rm: reconstMethodsVertLeft){
        wenoLevelsVertLeft.insert(rm.first);
    }

    mlrPtrVertLeft_->SelectWenoReconstLevel(wenoLevelsVertLeft, mlp);

    for (const auto& rm: reconstMethodsVertLeft){
        mlrPtrVertLeft_->ModifyReconstMethod(rm.first,rm.second);
    }

    mlrPtrVertLeft_->SeparateBoundaryLayer(mi);

    // Horizontal up method ============================================
    for (const auto& rm: reconstMethodsHoriUp){
        wenoLevelsHoriUp.insert(rm.first);
    }

    mlrPtrHoriUp_->SelectWenoReconstLevel(wenoLevelsHoriUp, mlp);

    for (const auto& rm: reconstMethodsHoriUp){
        mlrPtrHoriUp_->ModifyReconstMethod(rm.first,rm.second);
    }

    mlrPtrHoriUp_->SeparateBoundaryLayer(mi);

    // Horizontal down method ==========================================
    for (const auto& rm: reconstMethodsHoriDown){
        wenoLevelsHoriDown.insert(rm.first);
    }

    mlrPtrHoriDown_->SelectWenoReconstLevel(wenoLevelsHoriDown, mlp);

    for (const auto& rm: reconstMethodsHoriDown){
        mlrPtrHoriDown_->ModifyReconstMethod(rm.first,rm.second);
    }

    mlrPtrHoriDown_->SeparateBoundaryLayer(mi);
}

void NonSymDiffusion::diffusion::UpdateNonLinearWgts(const MeshInfo& mi){
    mlrPtrVertRight_->UpdateNonLinearWgts(mi,2);
    mlrPtrVertLeft_->UpdateNonLinearWgts(mi,2);
    mlrPtrHoriUp_->UpdateNonLinearWgts(mi,2);
    mlrPtrHoriDown_->UpdateNonLinearWgts(mi,2);
}

void NonSymDiffusion::diffusion::UpdateNonLinearWgtsAndDerivative(const MeshInfo& mi){
    mlrPtrVertRight_->UpdateNonLinearWgtsAndDerivs(mi);
    mlrPtrVertLeft_->UpdateNonLinearWgtsAndDerivs(mi);
    mlrPtrHoriUp_->UpdateNonLinearWgtsAndDerivs(mi);
    mlrPtrHoriDown_->UpdateNonLinearWgtsAndDerivs(mi);
}

void NonSymDiffusion::diffusion::GetInfo(const MeshInfo& mi){
    cout << "MLWENO reconstruction for Horizontal up side edges ..." << endl;
    mlrPtrHoriUp_->GetInfo();
    mlrPtrHoriUp_->PrintSmoothnessIndicator(mi);
    mlrPtrHoriUp_->PrintNonLinearWgts(mi);

    cout << "MLWENO reconstruction for Horizontal down side edges ..." << endl;
    mlrPtrHoriDown_->GetInfo();
    mlrPtrHoriDown_->PrintSmoothnessIndicator(mi);
    mlrPtrHoriDown_->PrintNonLinearWgts(mi);

    cout << "MLWENO reconstruction for Vertical right side edges ..." << endl;
    mlrPtrVertRight_->GetInfo();
    mlrPtrVertRight_->PrintSmoothnessIndicator(mi);
    mlrPtrVertRight_->PrintNonLinearWgts(mi);

    cout << "MLWENO reconstruction for Vertical left side edges ..." << endl;
    mlrPtrVertLeft_->GetInfo();
    mlrPtrVertLeft_->PrintSmoothnessIndicator(mi);
    mlrPtrVertLeft_->PrintNonLinearWgts(mi);

}
