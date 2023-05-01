#ifndef TRANSPORT_H_
#define TRANSPORT_H_

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

        unordered_set<int> boundaryCells;
        unordered_set<int> interiorCells;

        unordered_set<std::string> boundaryLevels;
        unordered_set<std::string> interiorLevels;

        // Create corresponding multilevel reconstruction with the given information.
        void CreateMLWENO(const MLWENO::MLWENOPrepare& mlp);

        // Update all flux and its derivatives (if implicit) on all the edges
        // collectively.
        void UpdateEdgeFlux(const MeshInfo& mi);

        void UpdateEdgeFluxDerivative(const MeshInfo& mi);

        // flag used to indicate if jacobian has been updated.
        int jacUpdate =0;

        // Return integrated flux corresponding to the given target cell
        double Flux(const MeshInfo& mi, const indice& global);

        // Compute derivatives of the integrated advection flux using in the jacobian.
        unordered_map<int, double> derivFlux(const MeshInfo& mi, 
                                             const indice& global);

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
         * Computation of flux around a given target cell.
         * Repeat calculation.
         * Should not be called directly.
         * Can be used for verification.
         * Collective update routine should be used instead.
         */
        double singleCellFlux_(const MeshInfo& mi, const indice& global, const double& t);
        unordered_map<int, double> singleCellDerivFlux_(const MeshInfo& mi, const indice& global, 
                                                        const double& t);

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

class diffusion {
    public:
        diffusion() {mlrPtr_ = new MLWENO::multiLevelReconstruction();};
        ~diffusion() {delete mlrPtr_;};

        /**
         * Diffusion defined on edges.
         * Different reconstruction methods used for horizontal and veritical edges.
         */
        unordered_map<std::string, vector<indice>> reconstMethodsVert;
        unordered_set<std::string> wenoLevelsVert;

        unordered_map<std::string, vector<indice>> reconstMethodsHori;
        unordered_set<std::string> wenoLevelsHori;

        unordered_set<int> boundaryCellsHori;
        unordered_set<int> interiorCellsHori;

        unordered_set<std::string> boundaryLevelsHori;
        unordered_set<std::string> interiorLevelsHori;

        unordered_set<int> boundaryCellsVert;
        unordered_set<int> interiorCellsVert;

        unordered_set<std::string> boundaryLevelsVert;
        unordered_set<std::string> interiorLevelsVert;

        std::string highestHoriLevel;
        std::string highestVertLevel;

        // Create corresponding multilevel reconstruction with the given information.
        void CreateMLWENO(const MLWENO::MLWENOPrepare& mlp);

    private:

        // Calculate diffusive flux.
        // And its corresponding derivatives.
        double flux_(const double * ru, int n);

        double dflux_(const double * ru, int n);

        // Interpolation positions
        const double alpha_ = 0.5;
        const double beta_ = 1.5;

        // Store calculated diffusive flux on the edge
        unordered_map<int, double> edgeHoriFlux_;
        unordered_map<int, double> edgeVertFlux_;

        /**
         * Pointer to a multilevel reconstruction class.
         */
        MLWENO::multiLevelReconstruction * mlrPtr_; 
};

class reaction {
    public:
        reaction() {mlrPtr_ = new MLWENO::multiLevelReconstruction();};
        ~reaction() {delete mlrPtr_;};

    private:
        MLWENO::multiLevelReconstruction * mlrPtr_;
};

class transport : public advection, public diffusion, public reaction {
    public:
        transport() {mlpPtr_ = new MLWENO::MLWENOPrepare();};
        ~transport() {delete mlpPtr_;};

        void AddLevel(const MeshInfo& mi, const int& stencilSizeX, 
                                          const int& stencilSizeY);

    private:

        /**
         * Holding all weno reconstruction in a pointer pointing to MLWENOPrepare object
         */
        MLWENO::MLWENOPrepare * mlpPtr_;
};

#endif
