#include "transportnew.h"

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
                                            const unordered_map<int, double>& duIn, 
                                            const unordered_map<int, double>& duOut){
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
    return flux_(uIn, fixed_, unitNormal, point, alphaLF);

    // Prescribe flux, Mimicing "Neumann" boundary condition
    //return 0.0;
}

unordered_map<int, double> boundaryCondition_(const double& uIn, const vertex& unitNOrmal,
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
double advection::Flux(const MeshInfo& mi, const indice& global, double t){

    double work = 0.0;

    //! Declare variable holding four corners of the given cell.
    vertexSet corner = extractCorners(mi, global); 

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

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

        double edgework = 0.0;

        //! Gauss quadrature rule.
        if (Interior_(mi,globalOut)){
            for (int g=0; g<gpe.size(); g++){
                vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

                double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,global);

                double uOut = mlrPtr_->EvaluateMLWENO(mi,mapped,globalOut);

                edgework += gwe[g] * flux_(uIn, uOut, unitNormal, mapped, uMax_) * len/2.0; 
            }
        }else{
             for (int g=0; g<gpe.size(); g++){
                vertex mapped = GaussMapPointsEdge({gpe[g]},edge);

                double uIn = mlrPtr_->EvaluateMLWENO(mi,mapped,global);

                edgework += gwe[g] * boundaryCondition_(uIn, unitNormal, mapped, uMax_) * len/2.0; 
            }
        }
        // Summing edge flux
        work += edgework;
    }

    work = work / NumIntegralFace(corner,{0,0}, {0.0,0.0}, 1.0, constFunc);

    return work;
}

/**
 * Compute derivative of advection flux
 */
unordered_map<int, double> advection::derivFlux(const MeshInfo& mi, const indice& global, 
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

// ========== Diffusion ===========================================

/**
 * Return diffusive flux
 */
double diffusion::flux_(const double * ru, int n){
    assert(n == 4);
    return ((ru[2]- ru[1])*beta_*beta_/(2*alpha_)-
            (ru[3]- ru[0])*alpha_*alpha_/(2*beta_))/
           (beta_*beta_-alpha_*alpha_);
}

double diffusion::dflux_(const double * dru, int n){
    assert(n == 4);
    return flux_(dru, n);
}

// ========== Reaction ============================================

// ========== Transport ===========================================
// Inherit from advection, diffusion and reaction class


