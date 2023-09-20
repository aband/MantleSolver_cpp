#ifndef MESHUSE_H_
#define MESHUSE_

#include "util.h" 
#include "input.h"
#include <petsc.h>

extern "C"{
#include "mesh.h"
#include "output.h"
}

class Initialize {
    public:
        //! A constructor
        /**!
         * Construct a Initialize class 
         * This Initialize class containing information regarding data management,
         * and mesh.
         */
        Initialize() {};

        //! A destructor
        /**!
         * Destruct a Initialize class.
         */
        ~Initialize() {VecDestroy(&globalu_); 
                       VecDestroy(&fullmesh_);
                       DMDestroy(&dmu_);
                       DMDestroy(&dmMesh_);};

        /**!
         * Function finish all the prepare work.
         */
        int Prepare();

        /**!
         * Initialize cell average values.
         */
        void CellAveragedInitialCondition(double (*func)(const valarray<double>& point, 
                                                         const vector<double>& param)); 

    private:
        // Global cell(element) size.
        int globalM_, globalN_;

        // Data management managing solution and mesh.
        DM dmu_, dmMesh_;

        // Vector storing full mesh.
        Vec fullmesh_;

        // Ghost layer size.
        int stencilWidthU_;
        int stencilWidthMesh_;

        double L_, H_;
        double xstart_, ystart_;

        // Test a single stencil for convergence
        int singleStencilTest_;
        double scale_;

        // Define type of mesh
        int meshtype_;

        Vec globalu_;

        friend class MeshUse;
};

class MeshUse {
    public:
        //! A constructor
        /**!
         * Construct a MeshUse class 
         * This MeshUse class containing meshInfo struct.
         * This MeshUse class manipulate members in meshInfo struct.
         * This MeshUse class is friend with MFEMUse and MLWENOUse class.
         * MeshInfo should not be called directly.
         */
        MeshUse(){};

        //! A destructor
        /**!
         * Destruct a MeshUse class.
         */
        ~MeshUse(){};

        /**!
         * Create mesh.
         */
        int CreateMeshInfo();

        /**
         * Finalize.
         */
        int Finalize();

    private:
        Vec localu_;

        MeshInfo mi_;

};
#endif
