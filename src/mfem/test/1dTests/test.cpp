#include <petsc.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include "integral.h"
#include "eutectic.h"
#include "input.h"
#include "util.h"

// MFEM parameter header file
#include "myFunc.h"

//#define COUPLED
#include "passemble.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "bndry.h"
#include "preconst.h"
#include "psolve.h"


extern "C"{
#include "mesh.h"
#include "output.h"
//#include "cgns_io.h"
}

int main(int argc, char ** argv){

    



    return 0;
}
