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

        //! A post work clean function
        /**!
         * Clean used dm and vec objects.
         * Should be called at the end of main function.
         */

        int clean();

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
        * And Mesh vector.
        */
        int CreateMesh(const int& M, const int& N,
                       double L, double H, 
                       double xstart, double ystart,
                       const int& stencilWidthMesh, 
                       const int& stencilWidthU,
                       const bool& physicsScale,
                       const int& meshType); 

        Vec globalmesh;  

        int PrintMesh();

       /**!
        * Assign Initial cell averaged condition.
        * Specificed for coupled system.
        * global vectors for dimensionless enthalpy and dimensionless Composition.
        */
        Vec globalCD, globalHD;
        Vec localCD, localHD;

        int InitCellAveVal(double (*funcHD)(const valarray<double>& point, const vector<double>& param),
                           double (*funcCD)(const valarray<double>& point, const vector<double>& param));

       // ================================================================================================
        


    private:
        /**!
         * Old file used in limited functions.
         */
        MeshParam mp_;

        /**!
         * WENO useage objects
         */
        MLWENO::MLWENOUse * mluse_;

        /**!
         * WENO preparation object.
         */
        MLWENO::MLWENOPrepare * mlpPtr_;

};

#endif
