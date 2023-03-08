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
            multiLevelReconstruction() {};

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
             * Call this function when a boundary layer size larger than 1 is needed.
             */
            void SpecifyBoundaryLayer(const MeshInfo& mi, const int& layerSize);

            //! Uniformally add reconstruction methods.
            void AddReconstMethod(vector<indice> brm) {baseReconstMethod_.push_back(brm); AddWgts_();}; 
            void AddReconstMethod(vector<vector<indice>> brms) 
                                 {for (int i=0; i<brms.size(); i++){
                                      AddReconstMethod(brms[i]);     
                                  }}; 

            //! Two ways of updateing non linear weights with updated smoothness indicators.
            void UpdateNonLinearWgts(const MeshInfo& mi);
            void UpdateTwoStageNonLinearWgts(const MeshInfo& mi);

            //! Routines used to check results and verification
            void GetInfo();
            void PrintSmoothnessIndicator(const MeshInfo& mi);
            void PrintNonLinearWgts(const MeshInfo& mi);

            void Clear();
        private:

            const double eps0_ = 0.01;

            map<std::string, singleLevelReconstruction *> reconstLevels_; //! reconstruction levels
            map<std::string, vector<indice> > reconstMethods_;            //! reconst methods

            vector< singleLevelReconstruction *> allLevels_;
            vector<vector<indice>> baseReconstMethod_; 
            vector< map<int, double> > linearWgts_;

            map<int, vector< map<int, double>>* > allLinearWgts_;
            map<int, vector< map<int, double> > > nonLinearWgts_;

            //! Bias usually set to be zero
            vector< map<int, int> > etaBias_;

            // Collective methods
            void ResetWgts_();
            void AddWgts_();

            void UpdateNonLinearWgts_(const MeshInfo& mi, indice start);
            void UpdateFirstStageNonLinearWgts_(const MeshInfo& mi, indice start);
            void UpdateTwoStageNonLinearWgts_(const MeshInfo& mi, indice start);

    };
}

#endif
