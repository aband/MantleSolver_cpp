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

/** !
 * Start with weakly coupled situation, where 
 * I approximate the Jacobian of the normal flux with diagonal matrix.
 */
namespace WEAK_COUPLED {

    enum bndryTypeTrans {"wall", "reflective", "absorb", "flux", "Dirichlet"}; 

    enum bndryTypeFlow {"essential", "natural"};
     
}

#endif
