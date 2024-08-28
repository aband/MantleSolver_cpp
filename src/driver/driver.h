#ifndef MESHUSE_H_
#define MESHUSE_H_

#include "util.h" 
#include "input.h"
#include <petsc.h>
#include "reconstMLWENO.h"

//#include "transport.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

enum transportType {adv, diff, adv_diff, adv_diff_react};

class Driver {
    public:
        //! A constructor
        /**!
         * Construct a driver class.
         * Driver class holding a pointer to meshInfo object.
         * Driver class will be used to interact with underlaying functions.
         */
        Driver() {};

        //! A destructor
        /**!
         * Destruct a Initialize class.
         */
        ~Driver() {};

        /**!
         * MeshInfo struct
         * Can be accessed from outside the class directly.
         */
        MeshInfo mi;

        /**!
         * Use WENO for reconstruction.
         * Allocate memory space for MLWENOPrepare class.
         */
        int UseWeno();

        /**!
         * Assign all possible stencil sizes to WENOPrepare class.
         */
        int AddLevel(const int& m, const int& n);

        /**!
         * Prepare for solving a transport problem 
         */
        int PrepareTransport(transportType type);

    private:

        /**!
         * Hold MLWENOPrepare pointer
         */
        MLWENO::MLWENOPrepare * mlpPtr_ = NULL;

        // Data management for mesh
        DM dmMesh_;       
        // Data management for solution
        DM dmu; 

};

#endif
