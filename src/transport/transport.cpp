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
 * Return diffusive flux
 */
double diffusion::flux(const double * ru, int n, double alpha, double beta){
    assert(n == 4);
    return ((ru[2]- ru[1])*beta*beta/(2*alpha)-
            (ru[3]- ru[0])*alpha*alpha/(2*beta))/
           (beta*beta-alpha*alpha);
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
double transport::advFlux(const MeshInfo& mi, const indice& global, double t){

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

            work += gwe[g] * advection::flux(uIn, uOut, unitNormal, mapped, 1.0) * len/2.0; 
        }

    }

    work = work / NumIntegralFace(corner,{0,0}, {0.0,0.0}, 1.0, constFunc);

    return work;
}

/**
 * Compute derivative of advection flux
 */
unordered_map<int, double> transport::derivAdvFlux(const MeshInfo& mi, const indice& global, 
                                                   const double& time){

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

    //! Compute area of target cell
    double area = NumIntegralFace(corner, {0,0}, {0.0,0.0}, 1.0, constFunc);
    //double area = mi.cellArea.at(FlatIndic(mi,global));

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
            unordered_map<int, double> derivIn = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,global);
            double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,global);

            unordered_map<int, double> derivOut;
            double uOut = advection::boundary;
            if (InsideBoundary_(mi, globalOut)){
                derivOut = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,global);
                uOut = mlrPtr_->EvaluateMLWENO(mi,mapped,global);
            }

            // Compute derivative of flux at a given gauss point
            unordered_map<int, double> derivflux = LaxFriedrichs::dflux(uIn, uOut, unitNormal,
                                                                        mapped, 1.0, 
                                                                        derivIn, derivOut);

            for (auto & derivf : derivflux){
                if (work.count(derivf.first) > 0){
                    work[derivf.first] += gwe[g]*derivf.second*len/2.0/area;
                } else {
                    work.insert(std::pair<int,double> (derivf.first, gwe[g]*derivf.second*len/2.0/area));
                }
            }

        } 
    }

    return work;
}

/**
 * Compute diffusion flux.
 * A compact version of function.
 */
double transport::edgeDiffFlux_(const MeshInfo& mi,
                                const indice& global, 
                                const vertexSet& edge,
                                const double& alpha,
                                const double& beta,
                                const double& scale){
    double work = 0.0;

    // Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Find Interpolation positions
    double len = length(edge);
    // Compute unit normal vector pointing outside.
    vertex unitNormal = UnitNormal(edge,len);

    for (int g = 0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        double ru[4];

        ru[0] = mlrPtr_->EvaluateMLWENO(mi,mapped-unitNormal*beta*scale, global);
        ru[1] = mlrPtr_->EvaluateMLWENO(mi,mapped-unitNormal*alpha*scale,global);
        ru[2] = mlrPtr_->EvaluateMLWENO(mi,mapped+unitNormal*alpha*scale,global);
        ru[3] = mlrPtr_->EvaluateMLWENO(mi,mapped+unitNormal*beta*scale, global);

        work += gwe[g] * diffusion::flux(ru, 4, alpha, beta) * len/2.0;

    }

    return work;
}

//! Cases where boundary values are fixed
double transport::edgeDiffFlux_(const MeshInfo& mi,
                                const indice& global, 
                                const vertexSet& edge,
                                const double& alpha,
                                const double& beta,
                                const double& scale,
                                const int * boundFix,
                                const int& n,
                                const double& boundaryValue){
    double work = 0.0;

    // Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Find Interpolation positions
    double len = length(edge);
    // Compute unit normal vector pointing outside.
    vertex unitNormal = UnitNormal(edge,len);

    for (int g = 0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        double ru[4];

        ru[0] = mlrPtr_->EvaluateMLWENO(mi,mapped-unitNormal*beta*scale, global);
        ru[1] = mlrPtr_->EvaluateMLWENO(mi,mapped-unitNormal*alpha*scale,global);
        ru[2] = mlrPtr_->EvaluateMLWENO(mi,mapped+unitNormal*alpha*scale,global);
        ru[3] = mlrPtr_->EvaluateMLWENO(mi,mapped+unitNormal*beta*scale, global);

        // Correct boundary values
        for (int i=0; i<n; i++){
            ru[boundFix[i]] = boundaryValue;
        }

        work += gwe[g] * diffusion::flux(ru, 4, alpha, beta) * len/2.0;

    }

    return work;
}

void transport::updateAllDiffFlux(const MeshInfo& mi){
    //! Calculate diffusive flux on edges.
    //! Calculated values are stored in diffusion::edgeFlux.
    // Clear previous data first.
    diffusion::edgeHoriFlux.clear();
    diffusion::edgeVertFlux.clear();

    double beta = diffusion::beta;
    double alpha = diffusion::alpha;

    int bF[2] = {0,0};

    const int n = 2;

    // In order to avoid repeating the same calculation
    // Calculate left and bottom edges for each cell first
    // Calculate right and top edges for the entire local part next 
    // Used predefined functions above to calculate left and bottom edges of each cell
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){
        indice globalCell {i,j};

        double scale = sqrt(mi.cellArea.at(FlatIndic(mi,globalCell)));

        // Extract four corners from mesh
        vertexSet corner = extractCorners(mi, globalCell);

        vertexSet horiEdge {corner.at(0), corner.at(1)};
        vertexSet vertEdge {corner.at(0), corner.at(3)};

        // Calculate diffusive flux
        // Calculate horizontal diffusive flux first
        if(j==0){
            // In this multi level calculation
            // The reconstruction of values exactly on the boundary will be 
            // using reconstruction methods one cell next to it
            bF[0] = 0;
            bF[1] = 1;
            indice shift {i,j+1};
            diffusion::edgeHoriFlux.insert(
                std::make_pair<int, double>(
                    FlatIndic(mi.MPIlocalCellSize[0], globalCell), 
                    edgeDiffFlux_(mi,shift,horiEdge,alpha,beta,scale,bF,n,diffusion::boundaryD)
                )
            );

        } else if (j==1 || j==mi.MPIlocalCellSize[1]-1){
            double tmpAlpha = alpha/beta;
            double tmpBeta = beta/beta;
            diffusion::edgeHoriFlux.insert(
                std::make_pair<int, double>(
                    FlatIndic(mi.MPIlocalCellSize[0], globalCell),
                    edgeDiffFlux_(mi,globalCell,horiEdge,tmpAlpha,tmpBeta,scale)
                )
            );
        } else {
            diffusion::edgeHoriFlux.insert(
                std::make_pair<int, double>(
                    FlatIndic(mi.MPIlocalCellSize[0], globalCell),
                    edgeDiffFlux_(mi,globalCell,horiEdge,alpha,beta,scale)
                )
            );
        }

        // Calculate vertical diffusive flux next
        if(i==0){
            // In this multi level calculation
            // The reconstruction of values exactly on the boundary will be 
            // using reconstruction methods one cell next to it
            bF[0] = 0;
            bF[1] = 1;
            indice shift {i+1,j};
            diffusion::edgeVertFlux.insert(
                std::make_pair<int, double>(
                    FlatIndic(mi.MPIlocalCellSize[0]+1, globalCell), 
                    edgeDiffFlux_(mi,shift,vertEdge,alpha,beta,scale,bF,n,diffusion::boundaryD)
                )
            );

        } else if (i==1 || i==mi.MPIlocalCellSize[0]-1){
            double tmpAlpha = alpha/beta;
            double tmpBeta = beta/beta;
            diffusion::edgeVertFlux.insert(
                std::make_pair<int, double>(
                    FlatIndic(mi.MPIlocalCellSize[0]+1, globalCell),
                    edgeDiffFlux_(mi,globalCell,vertEdge,tmpAlpha,tmpBeta,scale)
                )
            );
        } else {
            diffusion::edgeVertFlux.insert(
                std::make_pair<int, double>(
                    FlatIndic(mi.MPIlocalCellSize[0]+1, globalCell),
                    edgeDiffFlux_(mi,globalCell,vertEdge,alpha,beta,scale)
                )
            );
        }

    }}

    // Calculate remaining top and right edges on the boundary
    for(int i=0; i<mi.MPIlocalCellSize[0]; i++){
        indice globalCell {i,mi.MPIlocalCellSize[1]};

        indice shift {i,mi.MPIlocalCellSize[1]-1};

        bF[0] = 2;
        bF[1] = 3;

        double scale = sqrt(mi.cellArea.at(FlatIndic(mi,globalCell)));

        // Extract four corners from mesh
        vertexSet corner = extractCorners(mi, globalCell);

        vertexSet horiEdge {corner.at(0), corner.at(1)};

        diffusion::edgeHoriFlux.insert(
            std::make_pair<int, double>(
                FlatIndic(mi.MPIlocalCellSize[0], globalCell), 
                edgeDiffFlux_(mi,shift,horiEdge,alpha,beta,scale,bF,n,diffusion::boundaryD)
            )
        );
    }

    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
        indice globalCell {mi.MPIlocalCellSize[0],j};

        indice shift {mi.MPIlocalCellSize[0]-1,j};

        bF[0] = 2;
        bF[1] = 3;

        double scale = sqrt(mi.cellArea.at(FlatIndic(mi,globalCell)));

        // Extract four corners from mesh
        vertexSet corner = extractCorners(mi, globalCell);

        vertexSet vertEdge {corner.at(0), corner.at(3)};

        diffusion::edgeVertFlux.insert(
            std::make_pair<int, double>(
                FlatIndic(mi.MPIlocalCellSize[0]+1, globalCell), 
                edgeDiffFlux_(mi,shift,vertEdge,alpha,beta,scale,bF,n,diffusion::boundaryD)
            )
        );
    }

}

double transport::diffFlux(const MeshInfo& mi, const indice& global, const double& t){
    //! Combine diffusive flux on each edge of the target cell 

    double work = 0.0;

    //! Extracting diffusive flux on the boundary precalculated
    work += edgeHoriFlux.at(FlatIndic(mi.MPIlocalCellSize[0], global));
    work += edgeHoriFlux.at(FlatIndic(mi.MPIlocalCellSize[0], global[0], global[1]+1));

    work += edgeVertFlux.at(FlatIndic(mi.MPIlocalCellSize[0]+1, global));
    work += edgeVertFlux.at(FlatIndic(mi.MPIlocalCellSize[0]+1, global[0]+1, global[1]));

    work = work/mi.cellArea.at(FlatIndic(mi,global));

    return work;
}

/**
 * Compute derivative of diffusion flux using in the jacobian
 */
unordered_map<int, double> derivDiffFlux(const MeshInfo& mi, 
                                         const indice& global,
                                         const double& time){


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
