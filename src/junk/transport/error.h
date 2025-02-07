#ifndef ERROR_H_
#define ERROR_H_

#include "transport.h"
#include "input.h"

// Calculate L1 error at the given time
// if true solution is defined.

void L1Error(Vec *globalU,
             Vec *globalError,
             Vec *fullmesh,
             DM  dmu, 
             DM  dmMesh,
             const double& time,
             const MeshInfo& mi);

#endif
