#include "diffusion.h"

// ========== Diffusion ===========================================
// |      Containing functions regarding diffusion flux           |
// ================================================================

using namespace SymDiffusion;

/**
 * Return diffusive flux
 */
double diffusion::flux_(const double * ru, int n){
    assert(n == 4);
    return ((ru[2]- ru[1])*beta*beta/(2*alpha)-
            (ru[3]- ru[0])*alpha*alpha/(2*beta))/
           (beta*beta-alpha*alpha);
}

double diffusion::dflux_(const double * dru, int n){
    assert(n == 4);
    return flux_(dru, n);
}

/**
 * Check if a given target cell is inside the boundary
 * Return an integer that distinguishes different situations near the boundary.
 */
int diffusion::Interior_(const int& k, const int& size){

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

void diffusion::UpdateEdgeFlux(const MeshInfo& mi){

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

void diffusion::UpdateEdgeFluxDerivative(const MeshInfo& mi){
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

double diffusion::Flux(const MeshInfo& mi, const indice& global){

    double work = 0.0;

    double area = mi.cellArea.at(FlatIndic(mi,global));

    work += edgeHoriFlux_[FlatIndic(mi,global)]; 

    work += -1 * edgeHoriFlux_[FlatIndic(mi,global+mi.faceNormal[2])];

    work += -1 * edgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global+mi.faceNormal[1])];

    work += edgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global)];

    work /= area;
  
    return work;
}

unordered_map<int, double> diffusion::derivFlux(const MeshInfo& mi, const indice& global){

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

double diffusion::boundaryCondition_(double * ru, int n, const int& flag){

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
double diffusion::edgeFlux_(const MeshInfo& mi,
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

unordered_map<int, double> diffusion::derivEdgeFlux_(const MeshInfo& mi,
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


void diffusion::CreateMLWENO(const MLWENO::MLWENOPrepare& mlp, const MeshInfo& mi){
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

void diffusion::UpdateNonLinearWgts(const MeshInfo& mi){
    mlrPtrHori_->UpdateNonLinearWgts(mi,2);
    mlrPtrVert_->UpdateNonLinearWgts(mi,2);
}

void diffusion::GetInfo(const MeshInfo& mi){

    cout << "MLWENO reconstruction for Horizontal edges ..." << endl;
    mlrPtrHori_->GetInfo();
    mlrPtrHori_->PrintSmoothnessIndicator(mi);
    mlrPtrHori_->PrintNonLinearWgts(mi);

    cout << "MLWENO reconstruction for Vertical edges ..." << endl;
    mlrPtrVert_->GetInfo();
    mlrPtrVert_->PrintSmoothnessIndicator(mi);
    mlrPtrVert_->PrintNonLinearWgts(mi);
}

// New non symmetric weno reconstruction for diffusion flux



