#ifndef PHASECAL_H_
#define PHASECAL_H_

#include "phase.h"
#include <iostream>

double ComputePhase(const double& initXi, 
                    const double& hd,
                    const double& cbar,
                    EUTECTIC::evalPhase* pPtr);

#endif
