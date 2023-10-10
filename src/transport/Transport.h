#ifndef EDGEFLUX_H_
#define EDGEFLUX_H_

/**!
 * Combining advective and diffusive flux.
 * Presenting a total flux on the edge.
 */

#include "advectiveFlux.h"
#include "diffusiveFlux.h"
#include "fcns.h"

class Transport{

    public:
        /**!
         * Initialize data structure holding all flux on the edge.
         * flowTypes are "advection" "diffusion" and "adv-diff".
         */
        Transport(const MeshInfo& mi,
                  const std::array& flowType);
        ~Transport(){};

        /**!
         * Compute total flux on all the edges.
         */
        void getEdgeFlux(); 

    private:
        vector<double> edgeFlux_;
        bool isAdv = false;
        bool isDif = false;
};

#endif
