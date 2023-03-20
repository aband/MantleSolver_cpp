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
            double uOut = mlrPtr_->EvaluateMLWENO(mi,mapped,globalOut);

            work += gwe[g] * LaxFriedrichs::flux(uIn, uOut, unitNormal, mapped, 1.0) * len/2.0; 
        }
    }

    work = work / NumIntegralFace(corner,{0,0}, {0.0,0.0}, 1.0, constFunc);

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

void transport::Check(const MeshInfo& mi){
    mlrPtr_->GetInfo();

    mlrPtr_->PrintSmoothnessIndicator(mi);

    mlrPtr_->PrintNonLinearWgts(mi);
}
