#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "reconstruction.h"
#include "func.h"

class advection {
    public:
        advection() {};
        ~advection() {};

        unordered_map<std::string, vector<indice>> reconstMethods;
        unordered_set<std::string> wenoLevels;

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

        /**
         * Separate boundary layer
         */
        void SeparateBoundaryLayer(const MeshInfo& mi) {mlrPtr_->SeparateBoundaryLayer(mi);};

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
        double advFlux(const MeshInfo& mi, indice global, double t);

        /**
         * Check if there is anything wrong
         */
        void Check(const MeshInfo& mi);

    private:
        MLWENO::multiLevelReconstruction * mlrPtr_;
};

#endif
