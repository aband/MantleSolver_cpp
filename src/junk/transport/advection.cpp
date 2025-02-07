#include "advection.h"

// ========== Advection ===========================================

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
 * Return advection flux
 */ 
double advection::flux_(const double& uIn, const double& uOut, 
                        const vertex& unitNormal, const vertex& point){
    double work = 0.0;

    work = (funcX(point, uIn, 0) + funcX(point, uOut, 0))*unitNormal[0] + 
           (funcY(point, uIn, 0) + funcY(point, uOut, 0))*unitNormal[1]; 

    /**
     * Using local lax friedrichs stabilization without passing global factor.
     */
    double alphaLF = max(fabs(dfuncX(point, uIn, 0)*unitNormal[0] + dfuncY(point, uIn, 0)*unitNormal[1]),
                         fabs(dfuncX(point, uOut,0)*unitNormal[0] + dfuncY(point, uOut,0)*unitNormal[1]));

    work = 0.5 * (work - alphaLF*(uOut - uIn));

    return work;
}

double advection::flux_(const double& uIn, const double& uOut, 
                        const vertex& unitNormal, const vertex& point, 
                        const double& alphaLF){
    double work = 0.0;
        
    work = (funcX(point, uIn, 0) + funcX(point, uOut, 0))*unitNormal[0] + 
           (funcY(point, uIn, 0) + funcY(point, uOut, 0))*unitNormal[1]; 

    /** 
     * Passing global lax friedrichs stabilization factor into this function.
     */

    work = 0.5 * (work - alphaLF*(uOut - uIn));

    return work; 
}

unordered_map<int,double> advection::dflux_(const double& uIn, const double& uOut, const vertex& unitNormal, 
                                            const vertex& mapped, const double& alphaLF, 
                                            const derivative& duIn, 
                                            const derivative& duOut){
    unordered_map<int, double> work;

    double in = 0.5*(dfuncX(mapped, uIn, 0)*unitNormal[0] + 
                     dfuncY(mapped, uIn, 0)*unitNormal[1] + alphaLF);

    double out = 0.5*(dfuncX(mapped, uOut, 0)*unitNormal[0] + 
                      dfuncY(mapped, uOut, 0)*unitNormal[1] - alphaLF);

    //! Loop through derivative of outside cell
    for (auto& duin: duIn){
        // No need to check if key exists for the fact that work is now completely empty
        work.insert(std::pair<int, double>(duin.first, duin.second*in));
    }

    //! Check if it is outside the boundary
    if (duOut.empty()==0) {
        //! Loop through derivative of inside cell
        for (auto& duout: duOut){
            if (work.count(duout.first) > 0){
                // this key does exists
                work[duout.first] += duout.second*out;
            } else {
                // this key does not exists
                work.insert(std::pair<int,double>(duout.first, duout.second*out));
            }
        }
    }

    return work;
}

/** 
 * Check if a given cell is inside the boundary or not in the sense of 
 * multilevel reconstruction for advection.
 */
bool advection::Interior_(const MeshInfo& mi, const indice& target){

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (target[0] < 0 || target[0] > mi.MPIglobalCellSize[0] -1 ||
        target[1] < 0 || target[1] > mi.MPIglobalCellSize[1] -1 ){

        return false;
    } else {

        return true;
    }
}

/**
 * Define boundary condition for advection problem.
 * Allowing future definition of more complex boundary conditions.
 */
double advection::boundaryCondition_(const double& uIn, 
                                     const vertex& unitNormal, const vertex& point,
                                     const double& alphaLF){
    // Mimicing "Dirichlet" boundary condition
    //return flux_(uIn, fixed_, unitNormal, point, alphaLF);
    return flux_(uIn, fixed_, unitNormal, point, alphaLF);


    // Prescribe flux, Mimicing "Neumann" boundary condition
    //return 0.0;

}

unordered_map<int, double> advection::boundaryCondition_(const double& uIn, const vertex& unitNormal,
                                                         const vertex& mapped, const double& alphaLF,
                                                         const unordered_map<int,double>& duIn){
    unordered_map<int, double> tmp; 

    // Mimicing "Dirichlet" boundary condition
    // where derivative of uOut is zero.

    return dflux_(uIn, fixed_, unitNormal, mapped, uMax_, duIn, tmp);

    // Mimicing "Neumann" boundary condition.
    // where flux is prescribed on the boundary.
    // return an empty map is good enough.
    //return tmp;
}

// Integrated flux on edges with respect to the target cell
double advection::singleCellFlux(const MeshInfo& mi, const indice& global, const double& t){

    double work = 0.0;

    //! Declare variable holding four corners of the given cell.
    vertexSet corner = extractCorners(mi, global); 

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

        // Summing edge flux
        work += edgeFlux_(mi,global,globalOut,edge);
    }

    work = work / NumIntegralFace(corner,{0,0}, {0.0,0.0}, 1.0, constFunc);

    return work;
}

/**
 * Compute derivative of advection flux
 */
unordered_map<int, double> advection::singleCellDerivFlux(const MeshInfo& mi, const indice& global, 
                                                           const double& time){

    unordered_map<int, double> work;

    //! Declare variable holding four corners of the given cell.
    vertexSet corner = extractCorners(mi, global); 
    
    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

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
        if (Interior_(mi,globalOut)){
            for (int g=0; g<gpe.size(); g++){
                vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

                //! Compute derivative and value of multi level reconstruction
                unordered_map<int, double> derivIn = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,global);
                double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,global);

                unordered_map<int, double> derivOut = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,global);
                double uOut = mlrPtr_->EvaluateMLWENO(mi,mapped,global);

                // Compute derivative of flux at a given gauss point
                unordered_map<int, double> derivflux = dflux_(uIn, uOut, unitNormal,
                                                              mapped, uMax_, 
                                                              derivIn, derivOut);

                for (auto & derivf : derivflux){
                    if (work.count(derivf.first) > 0){
                        work[derivf.first] += gwe[g]*derivf.second*len/2.0/area;
                    } else {
                        work.insert(std::pair<int,double> (derivf.first, gwe[g]*derivf.second*len/2.0/area));
                    }
                }
            }
        } else {
            for (int g=0; g<gpe.size(); g++){
                vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

                //! Compute derivative and value of multi level reconstruction
                unordered_map<int, double> derivIn = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,global);
                double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,global);

                // Compute derivative of flux at a given gauss point
                unordered_map<int, double> derivflux = 
                    boundaryCondition_(uIn, unitNormal,mapped,uMax_, derivIn);

                // On flux confined boundary.
                // derivflux will be empty.
                // loop will not excute.
                for (auto & derivf : derivflux){
                    if (work.count(derivf.first) > 0){
                        work[derivf.first] += gwe[g]*derivf.second*len/2.0/area;
                    } else {
                        work.insert(std::pair<int,double> (derivf.first, gwe[g]*derivf.second*len/2.0/area));
                    }
                }
            }
        }
    }

    return work;
}

//! Collective flux update.
void advection::UpdateEdgeFlux(const MeshInfo& mi){

    // Udpate every left and bottom edge for each cell
    for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] + 1; j++){
    for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] + 1; i++){

        if (j<mi.MPIglobalCellSize[1] && i<mi.MPIglobalCellSize[0]){

        indice global {i,j};

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, global); 

        // Compute and restore Horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};

        //! Compute global indice of outside cell with respect to the inside cell.
        indice globalOut = global + mi.faceNormal[0];

        edgeHoriFlux_[FlatIndic(mi, global)] = edgeFlux_(mi, global, globalOut, hori);

        // Update top edge for the cells on the very top
        if (j == mi.MPIglobalCellSize[1]-1){
            hori = {corners.at(2), corners.at(3)};
            globalOut = global + mi.faceNormal[2];

            edgeHoriFlux_[FlatIndic(mi, globalOut)] = -1 * edgeFlux_(mi, global, globalOut, hori);
       }

        // Compute and restore Vertical flux
        vertexSet vert {corners.at(3), corners.at(0)};
        globalOut = global + mi.faceNormal[3];

        edgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,global)] = edgeFlux_(mi, global, globalOut, vert);

        // Update right edge for the cells on the very right of the local part
        if (i == mi.MPIglobalCellSize[0]-1){
            vert = {corners.at(1), corners.at(2)};
            globalOut = global + mi.faceNormal[1];
            edgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalOut)] = -1 * edgeFlux_(mi, global, globalOut, vert);
        }
 
        }
    }}

}

void advection::UpdateEdgeFluxDerivative(const MeshInfo& mi){

    // Udpate every left and bottom edge for each cell
    for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1]; j++){
    for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0]; i++){

        indice global {i,j};

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, global); 

        // Compute and restore Horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};

        //! Compute global indice of outside cell with respect to the inside cell.
        indice globalOut = global + mi.faceNormal[0];

        derivEdgeHoriFlux_[FlatIndic(mi, global)] = derivEdgeFlux_(mi, global, globalOut, hori);

        // Update top edge for the cells on the very top
        if (j == mi.MPIglobalCellSize[1]-1){
            hori = {corners.at(2), corners.at(3)};
            globalOut = global + mi.faceNormal[2];

            unordered_map<int, double> derivedgehori = derivEdgeFlux_(mi,global, globalOut, hori);
            for (auto& d: derivedgehori){
                d.second *= -1;
            }
            derivEdgeHoriFlux_[FlatIndic(mi, globalOut)] = derivedgehori;
       }

        // Compute and restore Vertical flux
        vertexSet vert {corners.at(3), corners.at(0)};
        globalOut = global + mi.faceNormal[3];

        derivEdgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,global)] = derivEdgeFlux_(mi, global, globalOut, vert);

        // Update right edge for the cells on the very right of the local part
        if (i == mi.MPIglobalCellSize[0]-1){
            vert = {corners.at(1), corners.at(2)};
            globalOut = global + mi.faceNormal[1];

            unordered_map<int,double> derivedgevert = derivEdgeFlux_(mi,global,globalOut,vert);
            for (auto& d: derivedgevert){
                d.second *= -1;
            }
            derivEdgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalOut)] = derivedgevert;
        }

    }}

    jacUpdate = 1;

}

void advection::UpdateEdgeFluxDerivativeFull(const MeshInfo& mi){

    // Udpate every left and bottom edge for each cell
    for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] ; j++){
    for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] ; i++){

        indice global {i,j};

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, global); 

        // Compute and restore Horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};

        //! Compute global indice of outside cell with respect to the inside cell.
        indice globalOut = global + mi.faceNormal[0];

        derivEdgeHoriFlux_[FlatIndic(mi, global)] = derivEdgeFluxFull_(mi, global, globalOut, hori);

        // Update top edge for the cells on the very top
        if (j == mi.MPIglobalCellSize[1]-1){
            hori = {corners.at(2), corners.at(3)};
            globalOut = global + mi.faceNormal[2];

            unordered_map<int, double> derivedgehori = derivEdgeFluxFull_(mi,global, globalOut, hori);
            for (auto& d: derivedgehori){
                d.second *= -1;
            }
            derivEdgeHoriFlux_[FlatIndic(mi, globalOut)] = derivedgehori;
       }

        // Compute and restore Vertical flux
        vertexSet vert {corners.at(3), corners.at(0)};
        globalOut = global + mi.faceNormal[3];

        derivEdgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,global)] = derivEdgeFluxFull_(mi, global, globalOut, vert);

        // Update right edge for the cells on the very right of the local part
        if (i == mi.MPIglobalCellSize[0]-1){
            vert = {corners.at(1), corners.at(2)};
            globalOut = global + mi.faceNormal[1];

            unordered_map<int,double> derivedgevert = derivEdgeFluxFull_(mi,global,globalOut,vert);
            for (auto& d: derivedgevert){
                d.second *= -1;
            }
            derivEdgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,globalOut)] = derivedgevert;
        }

    }}
} 

//! Calculate integrated flux defined on one given edge.
double advection::edgeFlux_(const MeshInfo& mi, 
                            const indice& globalIn,
                            const indice& globalOut,
                            const vertexSet& edge){

    double work = 0.0;

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Get edge lendth and unit vector normal to the given edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    int rank;   
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank); 

    if (Interior_(mi,globalOut)){

        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

            double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,globalIn);

            double uOut = mlrPtr_->EvaluateMLWENO(mi,mapped,globalOut);

            work += gwe[g] * flux_(uIn, uOut, unitNormal, mapped, uMax_) * len/2.0; 
        }

    } else {

        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

            double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,globalIn);

            work += gwe[g] * boundaryCondition_(uIn, unitNormal, mapped, uMax_) * len/2.0; 
        }

    }

    return work;
}

unordered_map<int,double> advection::derivEdgeFlux_(const MeshInfo& mi, 
                                                    const indice& globalIn,
                                                    const indice& globalOut,
                                                    const vertexSet& edge){

    unordered_map<int, double> work;

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);

    //! Compute unit normal vector pointing outside.
    vertex unitNormal = UnitNormal(edge,len);

    //! Loop through gauss points
    if (Interior_(mi,globalOut)){
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

            //! Compute derivative and value of multi level reconstruction
            unordered_map<int, double> derivIn = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,globalIn);
            double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,globalIn);

            unordered_map<int, double> derivOut = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,globalOut);
            double uOut = mlrPtr_->EvaluateMLWENO(mi,mapped,globalOut);

            // Compute derivative of flux at a given gauss point
            derivative derivflux = dflux_(uIn, uOut, unitNormal,
                                          mapped, uMax_, 
                                          derivIn, derivOut);

            double multi = gwe[g]*len/2.0;

            unordered_map_arithmetic(work,derivflux,std::plus<double>(),
                                          multi,std::multiplies<double>());
        }
    } else {
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

            //! Compute derivative and value of multi level reconstruction
            unordered_map<int, double> derivIn = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,globalIn);
            double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,globalIn);

            // Compute derivative of flux at a given gauss point
            derivative derivflux = boundaryCondition_(uIn, unitNormal,mapped,uMax_, derivIn);

            double multi = gwe[g]*len/2.0;

            unordered_map_arithmetic(work,derivflux,std::plus<double>(),
                                          multi,std::multiplies<double>());

        }
    }

    return work;
}

// Update derivative of flux fully including differentiation of non linear weights
unordered_map<int,double> advection::derivEdgeFluxFull_(const MeshInfo& mi, 
                                                        const indice& globalIn,
                                                        const indice& globalOut,
                                                        const vertexSet& edge){

    unordered_map<int, double> work;

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);

    //! Compute unit normal vector pointing outside.
    vertex unitNormal = UnitNormal(edge,len);

    derivative tmp;

    //! Loop through gauss points
    if (Interior_(mi,globalOut)){
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

            //! Compute derivative and value of multi level reconstruction
            unordered_map<int, double> derivIn = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,globalIn);
            tmp = mlrPtr_->EvaluateDerivMLWENOAdd(mi,mapped,globalIn);
            unordered_map_arithmetic(derivIn, tmp, std::plus<double>());
            tmp.clear();
            double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,globalIn);

            unordered_map<int, double> derivOut = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,globalOut);
            tmp = mlrPtr_->EvaluateDerivMLWENOAdd(mi,mapped,globalOut);
            unordered_map_arithmetic(derivOut, tmp, std::plus<double>());
            tmp.clear();
            double uOut = mlrPtr_->EvaluateMLWENO(mi,mapped,globalOut);

            // Compute derivative of flux at a given gauss point
            unordered_map<int, double> derivflux = dflux_(uIn, uOut, unitNormal,
                                                          mapped, uMax_, 
                                                          derivIn, derivOut);

            double modify = gwe[g]*len/2.0;

            unordered_map_arithmetic(work,derivflux,std::plus<double>(),modify,std::multiplies<double>());

        }
    } else {
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

            //! Compute derivative and value of multi level reconstruction
            unordered_map<int, double> derivIn = mlrPtr_->EvaluateDerivMLWENO(mi,mapped,globalIn);
            tmp = mlrPtr_->EvaluateDerivMLWENOAdd(mi,mapped,globalIn);
            unordered_map_arithmetic(derivIn, tmp, std::plus<double>());
            tmp.clear();
            double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,globalIn);

            // Compute derivative of flux at a given gauss point
            unordered_map<int, double> derivflux = 
                boundaryCondition_(uIn, unitNormal,mapped,uMax_, derivIn);

            double modify = gwe[g]*len/2.0;

            unordered_map_arithmetic(work,derivflux,std::plus<double>(),modify,std::multiplies<double>());

        }
    }

    return work;
}

//! Simply return values precalculated from collective update routine.
double advection::Flux(const MeshInfo& mi, const indice& global){

    double work = 0.0;

    double area = mi.cellArea.at(FlatIndic(mi,global));

    work += edgeHoriFlux_[FlatIndic(mi,global)]; 

    work += -1 * edgeHoriFlux_[FlatIndic(mi,global+mi.faceNormal[2])];

    work += -1 * edgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,global+mi.faceNormal[1])];

    work += edgeVertFlux_[FlatIndic(mi.MPIglobalCellSize[0]+1,global)];

    work /= area;

    return work;
}

unordered_map<int, double> advection::derivFlux(const MeshInfo& mi, const indice& global){

    derivative work;

    // Get scale variables relating to area of the target cell
    double area = mi.cellArea.at(FlatIndic(mi,global));

    double multi = -1.0/area;

    // bottom horizontal edge
    derivative tmp = derivEdgeHoriFlux_[FlatIndic(mi,global)];

    unordered_map_arithmetic(work, tmp, std::plus<double>(),
                                   area, std::divides<double>());

//    for (auto & derivf : tmp){
//        if (work.count(derivf.first) > 0){
//            work[derivf.first] += derivf.second/area;
//        } else {
//            work.insert(std::pair<int,double> (derivf.first, derivf.second/area));
//        }
//    }

    // top horizontal edge
    tmp = derivEdgeHoriFlux_[FlatIndic(mi,global+mi.faceNormal[2])];

    unordered_map_arithmetic(work, tmp, std::plus<double>(),
                                   multi, std::multiplies<double>());

//    for (auto & derivf : tmp){
//        if (work.count(derivf.first) > 0){
//            work[derivf.first] += -1*derivf.second/area;
//        } else {
//            work.insert(std::pair<int,double> (derivf.first, -1*derivf.second/area));
//        }
//    }

    // right vertical edge
    tmp = derivEdgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global+mi.faceNormal[1])];

    unordered_map_arithmetic(work, tmp, std::plus<double>(),
                                   multi, std::multiplies<double>());

//    for (auto & derivf : tmp){
//        if (work.count(derivf.first) > 0){
//            work[derivf.first] += -1*derivf.second/area;
//        } else {
//            work.insert(std::pair<int,double> (derivf.first, -1*derivf.second/area));
//        }
//    }

    // left vertical edge
    tmp = derivEdgeVertFlux_[FlatIndic(mi.MPIlocalCellSize[0]+1,global)];

    unordered_map_arithmetic(work, tmp, std::plus<double>(),
                                   area, std::divides<double>());

//    for (auto & derivf : tmp){
//        if (work.count(derivf.first) > 0){
//            work[derivf.first] += derivf.second/area;
//        } else {
//            work.insert(std::pair<int,double> (derivf.first, derivf.second/area));
//        }
//    }

    return work;
}

void advection::CreateMLWENO(const MLWENO::MLWENOPrepare& mlp, const MeshInfo& mi){

    assert(reconstMethods.empty() ==0);

    // Create weno levels from reconst method
    for (const auto& rm: reconstMethods){
        wenoLevels.insert(rm.first);
    }

    mlrPtr_->SelectWenoReconstLevel(wenoLevels,mlp); 

    for (const auto& rm: reconstMethods){
        mlrPtr_->ModifyReconstMethod(rm.first,rm.second);
    }

    mlrPtr_->SeparateBoundaryLayer(mi);
}

void advection::UpdateNonLinearWgts(const MeshInfo& mi){
    mlrPtr_->UpdateNonLinearWgts(mi,2);
}

void advection::UpdateNonLinearWgtsAndDerivative(const MeshInfo& mi){
    mlrPtr_->UpdateNonLinearWgtsAndDerivs(mi);
}

void advection::GetInfo(const MeshInfo& mi){
    mlrPtr_->GetInfo();
    mlrPtr_->PrintSmoothnessIndicator(mi);
    mlrPtr_->PrintNonLinearWgts(mi);
}
