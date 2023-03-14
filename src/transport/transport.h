#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "reconstruction.h"
#include "func.h"

class advection {
    public:
        advection() {};
        ~advection() {};

        void SelectReconstLevels(unordered_set<std::string> reconstLevels;);

        void Flux();

    private:
        unordered_set<std::string> reconstLevels_;

        /** 
         * Transport function and its derivative.
         * Returns function defined in the func.cpp file.
         */
        double Func_(double x, double y, double u, double t);
        double dFunc_(double x, double y, double y, double t);

};

class diffusion {
    public:
        diffusion() {};
        ~diffusion() {};

        void SelectReconstLevels(unordered_set<std::string> reconstLevels;);

        void Flux();

    private:
        unordered_set<std::string> reconstLevels_;
 
};

class reaction {
    public:
        reaction() {};
        ~reaction() {};

        void SelectReconstLevels(unordered_set<std::string> reconstLevels;);

    private:
        unordered_set<std::string> reconstLevels_;

}

class transport : public advection, public diffusion, public reaction {
    public:
        transport() {mlrPtr_ = new MLWENO::multiLevelReconstruction();};
        ~transport() {delete mlrPtr;};
  
        void AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY, vector<indice> brm);

    private:
        MLWENO::multiLevelReconstruction * mlrPtr_;
}

#endif
