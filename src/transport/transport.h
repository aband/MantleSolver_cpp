#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "reconstruction.h"
#include "func.h"

class advection {
    public:
        advection() {};
        ~advection() {};

        void Flux();

        unordered_map<std::string, vector<indice>> reconstMethods;
        unordered_set<std::string> wenoLevels;

    private:
        /** 
         * Transport function and its derivative.
         * Returns function defined in the func.cpp file.
         */
        double FuncX_(double x, double y, double u, double t);
        double dFuncX_(double x, double y, double y, double t);
        double FuncY_(double x, double y, double u, double t);
        double dFuncY_(double x, double y, double y, double t);
};

class diffusion {
    public:
        diffusion() {};
        ~diffusion() {};

        void Flux();

        unordered_map<std::string, vector<indice>> reconstMethods;
        unordered_set<std::string> wenoLevels;
};

class reaction {
    public:
        reaction() {};
        ~reaction() {};

        unordered_map<std::string, vector<indice>> reconstMethods;
        unordered_set<std::string> wenoLevels;
}

class transport : public advection, public diffusion, public reaction {
    public:
        transport() {mlrPtr_ = new MLWENO::multiLevelReconstruction();};
        ~transport() {delete mlrPtr;};

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
        unordered_set CreateWenoLevel(unordered_map<std::string, vector<indice>>& reconstMethods); 
    

    private:
        MLWENO::multiLevelReconstruction * mlrPtr_;
}

#endif
