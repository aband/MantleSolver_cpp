#ifndef PARAM_H_
#define PARAM_H_

#include "util.h"
#include "eutectic.h"
#include "myFunc.h"

using namespace EUTECTIC;

// Initial distribution of c_bar and h_bar
double ComputePorosity(const vertex& point, Phase * phase);

int PorosityOut(double xstart, double ystart, double L, double H, int seed, 
                Phase * phase);

// Transport boundary condition related

bool left_boundary(const indice& globalCell,
                   const MeshInfo& mi);

bool interior(const indice& globalCell, 
              const MeshInfo& mi);

bool edge(const indice& globalCell,
          const MeshInfo& mi);

bool corner(const indice& globalCell, 
            const MeshInfo& mi);

#endif
