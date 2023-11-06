#ifndef EDGEFLUX_H_
#define EDGEFLUX_H_

/**!
 * Combining advective and diffusive flux.
 * Presenting a total flux on the edge.
 */

#include "advectiveFlux.h"
#include "diffusiveFlux.h"
//#include "fcns.h"

enum flowType {advection, diffusion, adv_diff};

class EdgeFlux{

    public:
        /**!
         * Initialize data structure holding all flux on the edge.
         * flowTypes are "advection" "diffusion" and "adv-diff".
         */
        EdgeFlux(const MeshInfo& mi,
                 const flowType& fT);
        ~EdgeFlux(){};

        /**!
         * Compute total flux on all the edges.
         */
        void getEdgeFlux(const MeshInfo& mi,
                         const MLWENO::MLWENOUse& mluAdv,
                         const MLWENO::MLWENOUse& mluDif); 

    private:
        vector<double> edgeFlux_;
        bool isAdv = false;
        bool isDif = false;

        // advection stability coefficient
        double alpha_ = 1.0;

        // diffusion reconstruction scale
        double scale_ = 1.0;
};

#endif
