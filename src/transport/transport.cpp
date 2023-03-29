#include "transport.h"

//! Return function value defined in the separate function file
double advection::FuncX_(vertex x, double u, double t){
    return funcX(x, u, t);
}

double advection::dFuncX_(vertex x, double u, double t){
    return dfuncX(x, u, t);
}

double advection::FuncY_(vertex x, double u, double t){
    return funcY(x, u, t);
}

double advection::dFuncY_(vertex x, double u, double t){
    return dfuncY(x, u, t);
}

/**
 * Assign reconstruction method and weno levels to multi level reconstruction
 */
void transport::AssignReconstruction(const unordered_map<std::string, vector<indice>>& reconstMethods, const unordered_set<std::string>& wenoLevels){

    //! Modify reconstruction method first
    for (auto & pair: reconstMethods){
        mlrPtr_->ModifyReconstMethod(pair.first, pair.second);
    }

    mlrPtr_->SelectWenoReconstLevel(wenoLevels);
}

/**
 * Compute integrated advective flux.
 */
double transport::advFlux(const MeshInfo& mi, indice global, double t){

    double work = 0.0;

    //! Declare variable holding four corners of the given cell.
    vertexSet corner; 

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
     
    //! Retrieve local cell indice (including ghost vertex)
    indice ghostlayerShift {mi.vertexGhostLayerSize, mi.vertexGhostLayerSize};
    indice fullLocal = global - mi.MPIlocalCellStart + ghostlayerShift;

    //! Extract corners from mesh.
    for (auto & fcorner: mi.faceCorner){
        corner.push_back(mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0],fullLocal+fcorner)]);
    }

    //! Integral flux edge by edge.
    for (int pos = 0; pos < 4; pos ++){
        vertexSet edge;
        edge.push_back(corner[pos]);
        edge.push_back(corner[(pos+1)%4]);

        double len = length(edge);
        //! Compute unit normal vector pointing outside.
        vertex unitNormal = UnitNormal(edge,len);

        //! Compute global indice of outside cell with respect to the inside cell.
        indice globalOut = global + mi.faceNormal[pos];

        //! Gauss quadrature rule.
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

            double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,global);

            double uOut = 0.0;
            if (InsideBoundary_(mi,globalOut)){
                uOut = mlrPtr_->EvaluateMLWENO(mi,mapped,globalOut);
            }

            work += gwe[g] * LaxFriedrichs::flux(uIn, uOut, unitNormal, mapped, 1.0) * len/2.0; 
        }

    }

    work = work / NumIntegralFace(corner,{0,0}, {0.0,0.0}, 1.0, constFunc);

    return work;
}

/**
 * Compute derivative of advection flux
 */
const unordered_map<int, double>& derivAdvFlux(const MeshInfo& mi, const indice& global, 
                                               double time){

// loop through all reconstruction method
// if key exists add values
// if key not exists insert pair

    unordered_map<int, double> work;

    //! Declare variable holding four corners of the given cell.
    vertexSet corner; 

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
     
    //! Retrieve local cell indice (including ghost vertex)
    indice ghostlayerShift {mi.vertexGhostLayerSize, mi.vertexGhostLayerSize};
    indice fullLocal = global - mi.MPIlocalCellStart + ghostlayerShift;

    //! Extract corners from mesh.
    for (auto & fcorner: mi.faceCorner){
        corner.push_back(mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0],fullLocal+fcorner)]);
    }

    //! Loop through four edges of a given cell
    for (int pos=0; pos<4; pos++){

        vertexSet edge;
        edge.push_back(corner[pos]);
        edge.push_back(corner[(pos+1)%4]);

        double len = length(edge);
        //! Compute unit normal vector pointing outside.
        vertex unitNormal = UnitNormal(edge,len);

        //! Compute global indice of outside cell with respect to the inside cell.
        indice globalOut = global + mi.faceNormal[pos];

        //! Loop through gauss points
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

            //! Compute derivative and value of multi level reconstruction
            unordered_map<int, double> derivIn = mlrPtr_->EvaluateMLWENODeirv(mi,mapped,global);
            double uIn = mlrPtr->EvaluateMLWENO(mi,mapped,global);

            unordered_map<int, double> derivOut;
            double uOut = 0.0;
            if (InsideBoundary_(mi, globalOut)){
                derivOut = mlrPtr->EvaluateMLWENODeriv(mi,mapped,global);
                uOut = mlrPtr->EvaluateMLWENO(mi,mapped,global);
            }

           // ==========================????????????!!!!!!!!!!!!!!!!! 


        } 
    }

    return work;
}

// ==============================================================================================
void transport::AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY, vector<indice> brm){
    mlrPtr_->AddLevel(mi,stencilSizeX,stencilSizeY,brm);
}

void transport::AddReconstMethod(unordered_map<std::string, vector<indice>>& reconstMethods,
                                 std::string key, vector<indice> brm){
    //! Make sure the reconstruction key is new
    assert(reconstMethods.count(key) == 0);
    //! Add this new pair to reconstMethods
    reconstMethods.insert(std::pair<std::string, vector<indice>>(key, brm));
}

void transport::CreateWenoLevel(const unordered_map<std::string, vector<indice>>& reconstMethods, unordered_set<std::string>& wenoLevels){
    for (auto rm : reconstMethods){
        wenoLevels.insert(rm.first); 
    }
}

/** 
 * Check if a given cell is inside the boundary or not
 */
bool transport::InsideBoundary_(const MeshInfo& mi, const indice& target){

    if (target[0] < 0 || target[0] > mi.MPIglobalCellSize[0] -1 ||
        target[1] < 0 || target[1] > mi.MPIglobalCellSize[1] -1 ){
        return false;
    } else {
        return true;
    }

}

/**
 * Separate boundary layer for different sub problems
 */
void transport::SeparateAdvBoundaryLayer(const MeshInfo& mi){
    AssignReconstruction(advection::reconstMethods,
                         advection::wenoLevels);

    mlrPtr_->SeparateBoundaryLayer(mi);

    advection::boundaryCells = mlrPtr_->GetboundaryCells();
    advection::interiorCells = mlrPtr_->GetinteriorCells();

    advection::boundaryLevels = mlrPtr_->GetboundaryLevels();
    advection::interiorLevels = mlrPtr_->GetinteriorLevels();
}

void transport::SeparateDiffBoundaryLayer(const MeshInfo& mi){

    // Vertical edge reconstruction
    AssignReconstruction(diffusion::reconstMethodsVert,
                         diffusion::wenoLevelsVert);

    mlrPtr_->SeparateBoundaryLayer(mi);

    diffusion::boundaryCellsVert = mlrPtr_->GetboundaryCells();
    diffusion::interiorCellsVert = mlrPtr_->GetinteriorCells();

    diffusion::boundaryLevelsVert = mlrPtr_->GetboundaryLevels();
    diffusion::interiorLevelsVert = mlrPtr_->GetinteriorLevels();

    // Horizontal edge reconstruction
    AssignReconstruction(diffusion::reconstMethodsHori,
                         diffusion::wenoLevelsHori);

    mlrPtr_->SeparateBoundaryLayer(mi);

    diffusion::boundaryCellsVert = mlrPtr_->GetboundaryCells();
    diffusion::interiorCellsVert = mlrPtr_->GetinteriorCells();

    diffusion::boundaryLevelsVert = mlrPtr_->GetboundaryLevels();
    diffusion::interiorLevelsVert = mlrPtr_->GetinteriorLevels();
}

/**
 * Assign pre calculated boundary cells and levels to multi level reconstuction
 */
void transport::AssignBoundaryMethods(const unordered_set<int>& boundaryCells, 
                                      const unordered_set<int>& interiorCells,
                                      const unordered_set<std::string>& boundaryLevels, 
                                      const unordered_set<std::string>& interiorLevels){

    mlrPtr_->AssignboundaryCells(boundaryCells);
    mlrPtr_->AssigninteriorCells(interiorCells);

    mlrPtr_->AssignboundaryLevels(boundaryLevels);
    mlrPtr_->AssigninteriorLevels(interiorLevels);
}

void transport::Check(const MeshInfo& mi){
    mlrPtr_->GetInfo();

    mlrPtr_->PrintSmoothnessIndicator(mi);

    mlrPtr_->PrintNonLinearWgts(mi);
}
