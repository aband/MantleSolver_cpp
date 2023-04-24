#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "reconstruction.h"
#include "func.h"

class advection {
    public:
        advection() {};
        ~advection() {};

        // Pre computed reconstruction information
        // Will not be changed during computation
        unordered_map<std::string, vector<indice>> reconstMethods;
        unordered_set<std::string> wenoLevels;

        unordered_set<int> boundaryCells;
        unordered_set<int> interiorCells;

        unordered_set<std::string> boundaryLevels;
        unordered_set<std::string> interiorLevels;

        double flux(const double& uIn, const double& uOut, 
                    const vertex& unitNormal, const vertex& point,
                    const double& alphaLF) 
        {return LaxFriedrichs::flux(uIn, uOut, unitNormal, point, alphaLF);};

        unordered_map<int, double> dflux(const double& uIn, const double& uOut, const vertex& unitNormal, 
                                         const vertex& mapped, const double& alphaLF, 
                                         const unordered_map<int, double>& duIn, 
                                         const unordered_map<int, double>& duOut)
        {return LaxFriedrichs::dflux(uIn, uOut, unitNormal, mapped, alphaLF, duIn, duOut);};

    private:
        /** 
         * Transport function and its derivative.
         * Returns function defined in the func.cpp file.
         */
        double FuncX_(vertex x, double u, double t);
        double dFuncX_(vertex x, double u, double t);
        double FuncY_(vertex x, double u, double t);
        double dFuncY_(vertex x, double u, double t);
};

class diffusion {
    public:
        diffusion() {};
        ~diffusion() {};

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

        // Calculate diffusive flux
        double flux(const double * ru, int n, double alpha, double beta);

};

class reaction {
    public:
        reaction() {};
        ~reaction() {};

        unordered_map<std::string, vector<indice>> reconstMethods;
        unordered_set<std::string> wenoLevels;
};

class transport : public advection, public diffusion, public reaction {
    public:
        transport() {mlrPtr_ = new MLWENO::multiLevelReconstruction();};
        ~transport() {delete mlrPtr_;};

        /**
         * Compute possible weno levels in advance.
         * Add all useful weno levels to the class pointer.
         */
        void AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY, vector<indice> brm);

        /**
         * Define weno reconstruction method for different sub problems.
         */
        void AddReconstMethod(unordered_map<std::string, vector<indice>>& reconstMethods, 
                              std::string key, vector<indice> brm);

        /**
         * Create unordered set for weno levels with defined reconst method.  
         * Create this unordered set for the member function SelectWenoReconstLevel() from
         * multiLevelReconstruction pointer.
         */
        void CreateWenoLevel(const unordered_map<std::string, vector<indice>>& reconstMethods, unordered_set<std::string>& wenoLevels); 
  
        /**
         * Assign reconstruction method and weno levels to multi level reconstruction
         */
        void AssignReconstruction(const unordered_map<std::string, vector<indice>>& reconstMethods, const unordered_set<std::string>& wenoLevels);

        void AssignBoundaryMethods(const unordered_set<int>& boundaryCells, 
                                   const unordered_set<int>& interiorCells,
                                   const unordered_set<std::string>& boundaryLevels, 
                                   const unordered_set<std::string>& interiorLevels);

        /**
         * Separate boundary layer
         */
        void SeparateBoundaryLayer(const MeshInfo& mi) {mlrPtr_->SeparateBoundaryLayer(mi);};

        void SeparateAdvBoundaryLayer(const MeshInfo& mi);

        void SeparateDiffBoundaryLayer(const MeshInfo& mi);

        /**
         * Update non linear weights with given reconstruction methods and weno levels
         * Call update non linear weights from multi level reconstruction from under layer. 
         */
        void UpdateNonLinearWgts(const MeshInfo& mi, int stage) {mlrPtr_->UpdateNonLinearWgts(mi,stage);};

        /**
         * Compute advective flux.
         * Intended to write this function inside advection class.
         * Attempt failed.
         */
        double advFlux(const MeshInfo& mi, const indice& global, double t);

        /**
         * Compute derivative of advection flux using in the jacobian
         */
        unordered_map<int, double> derivAdvFlux(const MeshInfo& mi, 
                                                const indice& global,
                                                const double& time);

        /**
         * Compute diffusion flux.
         */
        double diffFlux(const MeshInfo& mi, const indice& global, const double& t);

        /**
         * Compute derivative of diffusion flux using in the jacobian
         */
        unordered_map<int, double> derivDiffFlux(const MeshInfo& mi, 
                                                 const indice& global,
                                                 const double& time);

        /**
         * Check if there is anything wrong
         */
        void Check(const MeshInfo& mi);

    private:

        /** 
         * Check if a given cell is inside the boundary or not
         */
        bool InsideBoundary_(const MeshInfo& mi, const indice& target);

        MLWENO::multiLevelReconstruction * mlrPtr_;
};

#endif
