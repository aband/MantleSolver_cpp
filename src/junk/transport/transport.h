#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "diffusion.h"
#include "advection.h"
#include "reaction.h"

using namespace NonSymDiffusion;

class transport : public advection, public diffusion, public reaction {
    public:
        transport() {mlpPtr_ = new MLWENO::MLWENOPrepare();};
        ~transport() {delete mlpPtr_;};

        /**
         * Add reconstruction levels to MLWENOPrepare class
         */
        void AddLevel(const MeshInfo& mi, const int& stencilSizeX, 
                                          const int& stencilSizeY);

        /**
         * Update smoothness indicator for all levels predefined.
         */
        void UpdateSmoothnessIndic(const MeshInfo& mi);

        void UpdateSmoothnessIndicAndDerivative(const MeshInfo& mi);

        /**
         * Assign reconstruction levels and methods to advection and diffusion
         */
        void AssignReconstMethod(unordered_map<std::string, vector<indice>>& reconstMethods,
                                 std::string key, const vector<indice>& brm);

        /**
         * Create MLWENO reconstruction for advection and diffusion class
         */
        void CreateMLWENO(const MeshInfo& mi);

    private:

        /**
         * Holding all weno reconstruction in a pointer pointing to MLWENOPrepare object
         */
        MLWENO::MLWENOPrepare * mlpPtr_;
};

#endif
