#ifndef ADVECTION_H_
#define ADVECTION_H_

#include "reconstruction.h"
#include "func.h"


/**
 * Incorporating MLWENOPrepare class.
 * Separate advection, diffusion and reaction multilevel reconstruction classes.
 * However, sharing the same MLWENOPrepare class.
 */

class advection {
    public:
        advection() {mlrPtr_ = new MLWENO::multiLevelReconstruction();};
        ~advection() {delete mlrPtr_;};

        // Pre computed reconstruction information
        // Will not be changed during computation
        unordered_map<std::string, vector<indice>> reconstMethods;
        unordered_set<std::string> wenoLevels;

        // Create corresponding multilevel reconstruction with the given information.
        void CreateMLWENO(const MLWENO::MLWENOPrepare& mlp, const MeshInfo& mi);

        // Update non linear weights
        void UpdateNonLinearWgts(const MeshInfo& mi);

        void UpdateNonLinearWgtsAndDerivative(const MeshInfo& mi);

        // Update all flux and its derivatives (if implicit) on all the edges
        // collectively.
        void UpdateEdgeFlux(const MeshInfo& mi);

        void UpdateEdgeFluxDerivative(const MeshInfo& mi);
        void UpdateEdgeFluxDerivativeFull(const MeshInfo& mi);

        // flag used to indicate if jacobian has been updated.
        int jacUpdate =0;

        // Return integrated flux corresponding to the given target cell
        double Flux(const MeshInfo& mi, const indice& global);

        // Compute derivatives of the integrated advection flux using in the jacobian.
        unordered_map<int, double> derivFlux(const MeshInfo& mi, 
                                             const indice& global);

        /**
         * Computation of flux around a given target cell.
         * Repeat calculation.
         * Should not be called directly.
         * Can be used for verification.
         * Collective update routine should be used instead.
         */
        double singleCellFlux(const MeshInfo& mi, const indice& global, const double& t);
        unordered_map<int, double> singleCellDerivFlux(const MeshInfo& mi, const indice& global, 
                                                       const double& t);
        // Print information
        void GetInfo(const MeshInfo& mi);

    private:
        /**
         * Define max value of u.
         * Used in global Lax-Friedrichs scheme as the stabilization constant.
         */
        const double uMax_ = 1.0;

        const double fixed_ = 0.0;

        /**
         * Define advection flux.
         * And its corresponding derivatives.
         * Several different flux scheme has been defined.
         * Flux will be automatically distinguished with the method of overloading.
         */
        //! Local Lax-Friedrichs flux scheme
        double flux_(const double& uIn, const double& uOut, 
                     const vertex& unitNormal, const vertex& point);

        //! Global Lax-Friedrichs flux scheme
        double flux_(const double& uIn, const double& uOut, 
                     const vertex& unitNormal, const vertex& point,
                     const double& alphaLF);

        //! Derivative of global Lax-Friedrich flux schme
        unordered_map<int, double> dflux_(const double& uIn, const double& uOut, const vertex& unitNormal, 
                                          const vertex& mapped, const double& alphaLF, 
                                          const unordered_map<int, double>& duIn, 
                                          const unordered_map<int, double>& duOut);

        /**
         * Calculate flux integral on one given edge.
         * Calculate derivative of the flux integral on one given edge at tha same time.
         */
        double edgeFlux_(const MeshInfo& mi, 
                         const indice& globalIn, 
                         const indice& globalOut, 
                         const vertexSet& edge);

        unordered_map<int,double> derivEdgeFlux_(const MeshInfo& mi, 
                                                 const indice& globalIn, 
                                                 const indice& globalOut, 
                                                 const vertexSet& edge);

        unordered_map<int,double> derivEdgeFluxFull_(const MeshInfo& mi, 
                                                     const indice& globalIn, 
                                                     const indice& globalOut, 
                                                     const vertexSet& edge);

        //! Store calculated advective flux and 
        //! corresponding derivatives on both horizontal and vertical edge
        //! The integer key in the following data structures are index of edges
        unordered_map<int, double> edgeHoriFlux_;
        unordered_map<int, double> edgeVertFlux_;

        unordered_map<int, unordered_map<int, double>> derivEdgeHoriFlux_;
        unordered_map<int, unordered_map<int, double>> derivEdgeVertFlux_;

        //! Determine whether the given target cell is inside the boudary or not.
        bool Interior_(const MeshInfo& mi, const indice& target);

        //! Deal with boundary condition
        //! Boundary condition returns boundary values for flux and derivatives of flux
        double boundaryCondition_(const double& uIn,  
                                  const vertex& unitNormal, const vertex& point,
                                  const double& alphaLF);

        unordered_map<int, double> boundaryCondition_(const double& uIn, const vertex& unitNormal, 
                                                      const vertex& mapped, const double& alphaLF, 
                                                      const unordered_map<int, double>& duIn);
        /** 
         * Transport function and its derivative.
         * Returns function defined in the func.cpp file.
         */
        double FuncX_(vertex x, double u, double t);
        double dFuncX_(vertex x, double u, double t);
        double FuncY_(vertex x, double u, double t);
        double dFuncY_(vertex x, double u, double t);

        /**
         * Pointer to a multilevel reconstruction class.
         */
        MLWENO::multiLevelReconstruction * mlrPtr_; 
};

#endif
