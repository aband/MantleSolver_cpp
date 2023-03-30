#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

//\slreconstruction.h slreconstruction.h
//! Weno reconstruction with a given stencil level

#include "polynomial.h"
#include <map>

namespace MLWENO{

    class singleLevelReconstruction {
        public:
            //! A constructor
            /*!
             * Default constructor for a single reconstruction level.
             */
            singleLevelReconstruction() {};
            //! A constructor
            /*!
             * Costume constructor with given x and y stencil size.
             */
            singleLevelReconstruction(int stencilSizeX, int stencilSizeY); 
           
            //! A destructor
            /*!
             * Clear interior, single level polynomial and smoothness indicators
             */
            ~singleLevelReconstruction();

            // ======================================================================
            void CreateStencilPolynomials(const MeshInfo& mi);

            bool CheckExist(const MeshInfo& mi, indice owner) const {return interior_.count(FlatIndic(mi,owner));};

            const int GetSizeX() const {return stencilSizeX_;};
            const int GetSizeY() const {return stencilSizeY_;};

            const double GetScale(int s) {return singleLevel_[s]->GetScale();};
            const double GetScale(const MeshInfo& mi, indice owner) {return singleLevel_[FlatIndic(mi,owner)]->GetScale();}; 

            //! Directly calculate smoothness indicator of a given stencil
            //! Should not be called directly for computational efficiency
            double CalculateSmoothnessIndic(const MeshInfo& mi, indice owner);

            /**
             * Directly extract pre-calculateed smoothness indicator.
             * Should always be the one to call when smoothness indicator is needed.
             */
            double GetSmoothnessIndic(const MeshInfo& mi, indice owner); 

            //! Update smoothness indicator for entire reconstruction level
            void UpdateSmoothnessIndic(const MeshInfo& mi);

            /**
             * Evaluate polynomial
             */
            double Evaluate(const MeshInfo& mi, indice owner, vertex point);

            // ======================================================================
            //! class members for checking and verification
            void CheckStencils() const {cout<< "Constructed "<< interior_.size() << " stencils with the size of " << stencilSizeX_ << " " << stencilSizeY_ << endl;};
            void CheckStencilPolynomials(const MeshInfo& mi, indice start);
            void PrintSmoothnessIndicator(const MeshInfo& mi);

        private:

            int stencilSizeX_ = -1; //!< Stencil size in x direction
            int stencilSizeY_ = -1; //!< Stencil size in y direction

            stencil <indice> stencilIndice_;           //!< Indices with given x and y sizes
            unordered_set<int> interior_;              //!< Numbering the created stencils
            map<int, double> smoothnessIndic_;         //!< Smoothness Indicators
            map<int, tensorProductPoly::stencilPolynomial*> singleLevel_; //!< Single level polynomials

            void IdentifyInteriorCell_(const MeshInfo& mi);

            const vertex ComputeStencilCenter_(const MeshInfo& mi, int flat);

            void ComputeStencilPolyn_(const MeshInfo& mi);
    };

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
             * Call the second function when a boundary layer size larger than 1 is needed.
             */
            void SeparateBoundaryLayer(const MeshInfo& mi, const int& layerSize);
            void SeparateBoundaryLayer(const MeshInfo& mi);

            //! Get a copy of calculated boundary and interior cells
            const unordered_set<int>& GetboundaryCells() const 
            {return boundaryCells_;};
            const unordered_set<int>& GetinteriorCells() const 
            {return interiorCells_;};

            //! Assing calculated boundary and interior cells
            void AssignboundaryCells (const unordered_set<int>& boundaryCells) 
            {boundaryCells_ = boundaryCells;};
            void AssigninteriorCells (const unordered_set<int>& interiorCells)
            {interiorCells_ = interiorCells;};

            //! Get a copy of calculated boundary and interior levels
            const unordered_set<std::string>& GetboundaryLevels() const 
            {return boundaryLevels_;};
            const unordered_set<std::string>& GetinteriorLevels() const 
            {return interiorLevels_;};

            //! Assing calculated boundary and interior levels
            void AssignboundaryLevels (const unordered_set<std::string>& boundaryLevels) 
            {boundaryLevels_ = boundaryLevels;};
            void AssigninteriorLevels (const unordered_set<std::string>& interiorLevels)
            {interiorLevels_ = interiorLevels;};

            /**
             * Update non linear weights for all levels.
             * Smoothness indicator will be updated every time when non linear weights updated.
             */
            void UpdateNonLinearWgts(const MeshInfo& mi, const int stage);

            /**
             * Update reconstruction methods.
             * Called when same level are used, but different reconstruction methods.
             * For example, advection and diffusion reconstructioin uses third order level 
             * but with different stencils.
             * When an existing level is added, update reconstruction method with new one.
             */
            void ModifyReconstMethod(std::string key, vector<indice> newRecosntMethod);
            void ModifyReconstMethod(int stencilSizeX, int stencilSizeY, vector<indice> newReconstMethod)
            {std::string key = '(' + std::to_string(stencilSizeX) + ',' + std::to_string(stencilSizeY) + ')';
             ModifyReconstMethod(key,newReconstMethod);};

            /**
             * Select Weno reconstruction levels to be used in the reconstruction.
             */
            void SelectWenoReconstLevel(unordered_set<std::string> keys);

            /**
             * Evaluation of given point with selected weno reconstruction method.
             */
            double EvaluateMLWENO(const MeshInfo& mi, vertex point, indice cell);

            /**
             * Derivative of the reconstruction of value with respect to the given point
             * with multi level weno method.
             */
            unordered_map<int, double> EvaluateDerivMLWENO(const MeshInfo& mi, const vertex& point, const indice& global);

            //! Routines used to check results and verification
            void GetInfo();
            void PrintBoundaryLayer(const MeshInfo& mi);
            void PrintSmoothnessIndicator(const MeshInfo& mi);
            void PrintNonLinearWgts(const MeshInfo& mi);

            /**
             * Clear everything.
             * Delete all pointers defined inside the map.
             * Can be called directly. 
             * May crush if destructor is called after Clear() is called.
             */
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

            unordered_map<int, unordered_map<std::string, unordered_map<int,double>>> nonLinearWgts_; //! Storing all non linear weights mapping to each cells

            /**
             * Separate different reconstruction levels for boundary and interior cells.
             * Immediately called after separating boundary layer.
             */
            void SeparateReconstMethods_();

            /**            
             * No need to define linear weights.
             * Simply average out the total number of levels; 
             */
            void UpdateOneStageNonLinearWgts_(const MeshInfo& mi, int flatGlobal, 
                                              const unordered_set<std::string>& levels);

            void UpdateTwoStageNonLinearWgts_(const MeshInfo& mi, int flatGlobal,
                                              const unordered_set<std::string>& levels); 

            //! Bias usually set to be zero
            vector< map<int, int> > etaBias_;

    };

}

#endif
