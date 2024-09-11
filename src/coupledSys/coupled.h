#ifndef COUPLED_H_
#define COUPLED_H_

// Transport
#include "mlwenouse.h"
#include "trans_setting.h"
#include "advectiveFlux.h"
#include "diffusiveFlux.h"

// Flow
#include "passemble.h"
#include "locmat.h"
#include "flow_setting.h"

// Phase
#include "eutectic.h"

/** !
 * Start with weakly coupled situation, where 
 * I approximate the Jacobian of the normal flux with diagonal matrix.
 */
namespace WEAK_COUPLED {

    enum bndryTypeTrans {"wall", "reflective", "absorb", "flux", "Dirichlet"}; 

    enum bndryTypeFlow {"essential", "natural"};

    /**!
     * Class holding a single transport.
     */
    class transport{
        public: 
            transport() {MLWENO::MLWENOUse * mluseAdv = new MLWENO::MLWENOUse();
				             MLWENO::MLWENOUse * mluseDif = new MLWENO::MLWENOUse();};
            ~transport(){delete mluseAdv;
                         delete mluseDif};

            bool advFlag = false;
            bool difFlag = false;

            MLWENO::MLWENOUse * mluseAdv = NULL;
            MLWENO::MLWENOUse * mluseDif = NULL;

    }

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

    };



}

#endif
