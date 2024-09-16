#ifndef COUPLED_H_
#define COUPLED_H_

// Transport
#include "mlwenouse.h"
#include "advectiveFlux.h"
#include "diffusiveFlux.h"

// Phase
#include "eutectic.h"

using namespace EUTECTIC; 

/**!
 * Transport location functions
 * A set of functions define location of a given cell.
 */

bool interior(const indice& globalCell, 
              const MeshInfo& mi);

bool edge(const indice& globalCell,
          const MeshInfo& mi);

bool corner(const indice& globalCell, 
            const MeshInfo& mi);

/** !
 * Start with weakly coupled situation, where 
 * I approximate the Jacobian of the normal flux with diagonal matrix.
 */
namespace WEAK_COUPLED {

    /**!
     * Class holding a single transport.
     */
    class transport{
        public: 
            transport() {MLWENO::MLWENOUse * mluseAdv = new MLWENO::MLWENOUse();
                         MLWENO::MLWENOUse * mluseDif = new MLWENO::MLWENOUse();};
            ~transport(){delete mluseAdv;
                         delete mluseDif;};

            bool advFlag = false;
            bool difFlag = false;

            MLWENO::MLWENOUse * mluseAdv = NULL;
            MLWENO::MLWENOUse * mluseDif = NULL;

            std::string location(const MeshInfo& mi, const indice& globalCell);

            double reconstVal(const vertex& mapped, 
                              const indice& globalCell,
                              const MeshInfo& mi) 
            {return mluseAdv->Evaluate(mapped, globalCell, mi, location(mi, globalCell));};

    };

    /** !
     * Class holding all mlwenouse objects that used in the simulation of 
     * coupled transport.
     */
    class coupledTrans{
        public:
            coupledTrans() {transport * transCompon = new transport();
                            transport * transEnthal = new transport();};
            ~coupledTrans() {delete transCompon;
                             delete transEnthal;};

            transport * transCompon = NULL;
            transport * transEnthal = NULL;


            void reconstVal(const MeshInfo& mi, 
                            const indice& globalCell,
                            const vertex& mapped) 
            {HD = transEnthal->reconstVal(mapped, globalCell, mi);
             CD = transCompon->reconstVal(mapped, globalCell, mi);};

            // Current reconstructed values of nondimensional enthalpy and component
            double HD;
            double CD;
    };
}

#endif
