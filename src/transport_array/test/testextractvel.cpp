#include "petsc.h"
#include "input.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

#include "tensorstencilpoly.h"
#include "reconstruction.h"
#include <chrono>

#include "error.h"

#include "advectiveflux.h"

#include "rk.h"

#include "extractvel.h"

int main(int argc, char ** argv){

    vector<vertex> testfield;

    // A 3*3 cells test case
    constVelField(testfield, 3, 3);

    return 1;
}
