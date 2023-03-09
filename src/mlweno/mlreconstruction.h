#ifndef MLRECONSTRUCTION_H_
#define MLRECONSTRUCTION_H_

#include "slreconstruction.h"

namespace MLWENO{
    class multiLevelReconstruction {
        public:
            //! A constructor
            /**!
             * Construt multi-level weno reconstruction by specifying each single level
             */
            //multiLevelReconstruction() {};

            multiLevelReconstruction(const MeshInfo& mi, int stencilSizeX, int stencilSizeY, 
                                     vector<indice> brm)
            {AddLevel(mi,stencilSizeX,stencilSizeY,brm);};

            multiLevelReconstruction(const MeshInfo& mi, int* stencilSize, vector<indice> brm) 
            {AddLevel(mi,stencilSize,brm);};

            multiLevelReconstruction(const MeshInfo& mi, vector<int*> stencilSizes,
                                     vector<vector<indice>> brms) 
            {AddLevel(mi,stencilSizes,brms);};

            //! A destructor
            /**
             * Construct multi-level weno reconstruction by specifying each single level
             */
            ~multiLevelReconstruction() {Clear();};

            /**
             * Add a single weno reconstruction level to the weno reconstruction.
             * Create a reference key name for specific level.
             * Create map from the reference key to the single level reconstruction.
             */
            void AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY, vector<indice> brm);
            void AddLevel(const MeshInfo& mi, int* stencilSize, vector<indice> brm) {AddLevel(mi,stencilSize[0],stencilSize[1],brm);};
            void AddLevel(const MeshInfo& mi, vector<int*> stencilSizes, vector<vector<indice>> brms)
            {for (int i=0; i<stencilSizes.size(); i++){
                 AddLevel(mi,stencilSizes[i],brms[i]);}};

            /**
             * Specify boundary layers (cells near boundary that need additional reconstruciton level then interior cells)
             * Default boundary layer size is set to be 1.
             * Call the second function when a boundary layer size larger than 1 is needed.
             */
            void SeparateBoundaryLayer(const MeshInfo& mi, const int& layerSize);
            void SeparateBoundaryLayer(const MeshInfo& mi);

            /**
             * Update non linear weights for all levels.
             * Smoothness indicator will be updated every time when non linear weights updated.
             */
            void UpdateOneStageNonLinearWgts(const MeshInfo& mi);
            void UpdateTwoStageNonLinearWgts(const MeshInfo& mi);

            //! Routines used to check results and verification
            void GetInfo();
            void PrintBoundaryLayer(const MeshInfo& mi);
            void PrintSmoothnessIndicator(const MeshInfo& mi);
            void PrintNonLinearWgts(const MeshInfo& mi);

            void Clear();
        private:

            const double eps0_ = 0.01;

            std::string lowestLevel_;
            std::string highestLevel_;

            map<std::string, singleLevelReconstruction *> reconstLevels_; //! reconstruction levels
            map<std::string, vector<indice> > reconstMethods_;            //! reconst methods

            unordered_set<int> interiorCells_; //! Cells inside the computational domain.
            unordered_set<int> boundaryCells_; //! Cells on the boundary layer requiring additional resolution.

            unordered_set<std::string> wenoLevels_;     //! Storing all reconstruction levels defined previously.
            unordered_set<std::string> interiorLevels_; //! Reconstruction levels for interior cells
            unordered_set<std::string> boundaryLevels_; //! Reconstruction levels for boundary cells

            int totalLevels_ = 0; //! Accumulate all number of stencils.

            /**
             * Separate different reconstruction levels for boundary and interior cells.
             * Immediately called after separating boundary layer.
             */
            void SeparateReconstMethods_();

            /**            
             * No need to define linear weights.
             * Simply average out the total number of levels; 
             */
            unordered_map<int, unordered_map<std::string, unordered_map<int,double>>> nonLinearWgts; //! 

            //! Bias usually set to be zero
            vector< map<int, int> > etaBias_;

    };
}

#endif
