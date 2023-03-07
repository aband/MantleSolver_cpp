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
            multiLevelReconstruction(const MeshInfo& mi, int stencilSizeX, int stencilSizeY){AddLevel(mi,stencilSizeX,stencilSizeY);};
            multiLevelReconstruction(const MeshInfo& mi, int* stencilSize) {AddLevel(mi,stencilSize);};
            multiLevelReconstruction(const MeshInfo& mi, vector<int*> stencilSizes) 
            {AddLevel(mi,stencilSizes);};

            ~multiLevelReconstruction() {Clear();};

            //! A destructor
            /**!
             * Construt multi-level weno reconstruction by specifying each single level
             */
            //! Uniformally add reconstruction level.
            void AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY);
            void AddLevel(const MeshInfo& mi, int* stencilSize) {AddLevel(mi,stencilSize[0],stencilSize[1]);};
            void AddLevel(const MeshInfo& mi, vector<int*> stencilSizes)
            {for (int i=0; i<stencilSizes.size(); i++){
                 AddLevel(mi,stencilSizes[i]);}};

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

            vector< singleLevelReconstruction *> allLevels_;
            vector<vector<indice>> baseReconstMethod_; 
            vector< map<int, double> > linearWgts_;

            map <int, vector<vector<indice>>* > reconstMethods_;
            map <int, vector< map<int, double>>* > allLinearWgts_;
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
