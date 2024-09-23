#ifndef DRIVER_H_
#define DRIVER_H_

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

#define COUPLED
#include "passemble.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "bndry.h"
#include "preconst.h"
#include "psolve.h"

// MLWENO parameter header file
#include "mlwenouse.h"
#include "coupled.h"

extern "C"{
#include "mesh.h"
#include "output.h"
//#include "cgns_io.h"
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
        ~Driver() {delete myPhase;};

        /**!
         * MeshInfo struct
         * Can be accessed from outside the class directly.
         */
        MeshInfo mi;

       /**!
         * Data management objects for mesh and solution.
         * showing up in compuation process.
         */ 
        DM dmMesh;       
        DM dmu; 

       /**!
		  * Initialize phase package
		  */
       Phase * myPhase;

       int CreatePhase();

       /**!
		  * Create Data management objects.
		  */
        int CreateDMs(const int& M, const int& N,
                      double L, double H, 
							 double xstart, double ystart,
                      const int& stencilWidthMesh, 
							 const int& stencilWidthU,
							 const bool& physicsScale); 

        int PrintMesh();

       /**!
		  * Create Mesh vector
		  */

        Vec 
        int CreateMesh(); 

    private:
        /**!
         * Old file used in limited functions.
         */
        MeshParam mp_;

};

#endif
