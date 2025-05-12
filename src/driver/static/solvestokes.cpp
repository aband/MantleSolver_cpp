#include "driver.h"

int Driver::solveStokes(int maxIter, double tolUzawa){

    ParallelMatrixAssemble_case(allwgts, lphi);

    int nelem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    CreateLinearSys(reducedStokes_, nelem);
    CreateLinearSys(reducedDarcy_, nelem);

    Uzawa(reducedStokes_,nelem); 

    return 1;
}
