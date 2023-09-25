#ifndef MESHUSE_H_
#define MESHUSE_

#include "util.h" 
#include "input.h"
#include <petsc.h>

extern "C"{
#include "mesh.h"
#include "output.h"
}

class Driver {
    public:
        //! A constructor
        /**!
         * Construct a driver class.
         * Driver class holding a pointer to meshInfo object.
         * Driver class will be used to interact with underlaying functions.
         */
        Driver(MeshInfo * mi) {mi = &mi_;};

        //! A destructor
        /**!
         * Destruct a Initialize class.
         */
        ~Driver() {};

        /**!
         * Prepare for 
         */
        int PrepareTransport();

    private:
        // MeshInfo struct
        MeshInfo mi_;

        // Transport pointer

};

#endif
