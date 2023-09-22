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
         * Construct a Initialize class 
         * This Initialize class containing information regarding data management,
         * and mesh.
         */
        Driver() {};

        //! A destructor
        /**!
         * Destruct a Initialize class.
         */
        ~Driver();

        /**!
         * Function finish all the prepare work.
         */
        int Prepare();

        /**!
         * Initialize cell average values.
         */
        void CellAveragedInit(double (*func)(const valarray<double>& point, 
                                             const vector<double>& param)); 

        /**!
         * Create MeshInfo object.
         */
        int CreateMeshInfo();

    private:
        // Storing global cell size.
        int globalM_, globalN_;

        // Storing physical domain size.
        double L_, H_;

        // Data management managing solution and mesh.
        DM dmu_, dmMesh_;

        // Global vector storing full mesh.
        Vec fullmesh_;

        // Global and local vector storing cell averaged solution.
        Vec globalu_, localu_;

        // MeshInfo struct
        MeshInfo mi_;
};

#endif
