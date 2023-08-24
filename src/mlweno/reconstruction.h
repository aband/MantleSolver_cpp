#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

#include "polynomial.h"
#include <map>

namespace MLWENO{

    //! Single level reconstruction class
    /**
     * 1. Create non-overlap stencils with given stencil size.
     * 2. Compute stencil polynomials with created stencils. 
     * 3. Calculate smoothness indicators with stencil polynomials.
     */
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

            void CreateStencilPolynomials(const MeshInfo& mi);

            bool CheckExist(const MeshInfo& mi, indice owner) const {return interior_.count(FlatIndic(mi,owner));};

            const int GetSizeX() const {return stencilSizeX_;};
            const int GetSizeY() const {return stencilSizeY_;};

            const double GetScale(int s) {return singleLevel_[s]->GetScale();};
            const double GetScale(const MeshInfo& mi, const indice& owner) {return singleLevel_[FlatIndic(mi,owner)]->GetScale();}; 
            //! Directly calculate smoothness indicator of a given stencil
            //! Should not be called directly for computational efficiency
            double CalculateSmoothnessIndic(const MeshInfo& mi, indice owner);

            /**
             * Directly extract pre-calculateed smoothness indicator.
             * Should always be the one to call when smoothness indicator is needed.
             */
            double GetSmoothnessIndic(const MeshInfo& mi, indice owner); 

            unordered_map<int, double> GetSmoothnessIndicDeriv(const MeshInfo& mi, const indice& owner);

            //! Update smoothness indicator for entire reconstruction level
            void UpdateSmoothnessIndic(const MeshInfo& mi);

            void UpdateDerivSmoothnessIndic(const MeshInfo& mi);
            /**
             * Evaluate polynomial
             */
            double Evaluate(const MeshInfo& mi, const indice& owner, const vertex& point);

            double Evaluate(const MeshInfo& mi, const indice& owner, const vertex& point, const int& local);

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
            map<int, double> smoothnessIndic_;         //!< Smoothness indicators
            unordered_map<int, derivative> smoothnessIndicDeriv_;//!< Derivative of smoothness indicators
            map<int, tensorProductPoly::stencilPolynomial*> singleLevel_;       //!< Single level polynomials

            void IdentifyInteriorCell_(const MeshInfo& mi);

            const vertex ComputeStencilCenter_(const MeshInfo& mi, int flat);

            void ComputeStencilPolyn_(const MeshInfo& mi);
    };

    /**
     * A class preparing for all possible multi-level weno reconstruciton using
     * different single-level reconstruction.
     * This class holding all single level reconstuctions.
     * Let MLWENO "steal" single level reconstructions from it.
     * Multi-level reconstruction class is declared as a friend of this class.
     */
    class MLWENOPrepare {
        public:
            MLWENOPrepare() {};

            ~MLWENOPrepare();

            /**
             * Add single levels to the private allLevels_ member.
             * Reconstruction method not required.
             * Reconstruction method will be defined later in multiLevelReconstruction.
             */
            void AddLevel(const MeshInfo& mi, const int& stencilSizeX, 
                                              const int& stencilSizeY);

            /**
             * Update smoothness indicator for all single level stencil polynomials.
             * Should be called everytime non linear weights are being calculated.
             */
            void UpdateSmoothnessIndic(const MeshInfo& mi);

            /**
             * Update derivatives of smoothness indicator for all single level 
             * stencil polynomials.
             * Should be called before using fully implicit method.
             */
            void UpdateDerivSmoothnessIndic(const MeshInfo& mi);

            /**
             * Print information of all levels created.
             */
            void PrintInfo();

        private:
            friend class multiLevelReconstruction;

            std::unordered_map<std::string, singleLevelReconstruction *> allLevels_;
    };

    class multiLevelReconstruction {
        public:
            //! A constructor
            /**!
             * Construt multi-level weno reconstruction by specifying each single level
             */
            multiLevelReconstruction() {};

            /**
             * Be careful with the parametrized constructor function.
             * It is better to "steal" from MLWENOPrepare, instead of calculate itself.
             */
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

            void SeparateBoundaryLayer(const MeshInfo& mi, const int& layerSize, 
                                       const unordered_set<std::string>& additionalLevels); 

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
             * Update non linear weights default with two stage method
             */
            void UpdateNonLinearWgts(const MeshInfo& mi);

            void UpdateNonLinearWgtsAndDerivs(const MeshInfo& mi);

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
            void SelectWenoReconstLevel(const unordered_set<std::string>& keys);

            //! A completelly different function despite sharing the same name with
            //! the previous one. "Stealing" single level reconstruction from 
            //! class MLWENOPrepare.
            void SelectWenoReconstLevel(const unordered_set<std::string>& keys,
                                        const MLWENOPrepare& mlpPtr);

            /**
             * Evaluation of given point with selected weno reconstruction method.
             */
            double EvaluateMLWENO(const MeshInfo& mi, vertex point, indice cell) const;

            /**
             * Derivative of the reconstruction of value with respect to the given point
             * with multi level weno method.
             */
            unordered_map<int, double> EvaluateDerivMLWENO(const MeshInfo& mi, 
                                                           const vertex& point, 
                                                           const indice& global) const;

            unordered_map<int, double> EvaluateDerivMLWENOAdd(const MeshInfo& mi,
                                                              const vertex& point, 
                                                              const indice& global) const;
            /**
             * Get scale of the selected single level reconstruction
             * with respect to the selected cell.
             */
            const double GetScale(std::string level, int s)
            {return reconstLevels_[level]->GetScale(s);};
            const double GetScale(std::string level, const MeshInfo& mi, const indice& owner)
            {return reconstLevels_[level]->GetScale(mi,owner);};

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

            //! Storing all non linear weights mapping to each cells
            unordered_map<int, unordered_map<std::string, unordered_map<int,double>>> nonLinearWgts_; 

            //! Storing all derivatives of non liear weights mapping to each cells
            // The most outside map creates map from global index of cells to a set of non linear weights
            // The inner map from levels (std::string) to a set of non linear weights (within the same level)
            // The most inside map connects derivatives (a map mapping with cells in the chosen stencil).
            unordered_map<int, unordered_map<std::string, unordered_map<int,unordered_map<int,double>>>> derivNonLinearWgts_;

            /**
             * Separate different reconstruction levels for boundary and interior cells.
             * Immediately called after separating boundary layer.
             */
            void SeparateReconstMethods_();

            void SeparateReconstMethods_(const unordered_set<std::string>& additionalLevels);

            /**            
             * No need to define linear weights.
             * Simply average out the total number of levels; 
             */
            void UpdateOneStageNonLinearWgts_(const MeshInfo& mi, int flatGlobal, 
                                              const unordered_set<std::string>& levels);

            void UpdateTwoStageNonLinearWgts_(const MeshInfo& mi, int flatGlobal,
                                              const unordered_set<std::string>& levels); 

            void UpdateNonLinearWgts_(const MeshInfo& mi, int flatGlobal,
                                      const unordered_set<std::string>& levels,
                                      const unordered_map<std::string, int>& powerShift);

            /**
             * Update non linear weights and its 
             * corresponding derivatives at the same time.
             */
            void UpdateNonLinearWgtsAndDerivs_(const MeshInfo& mi, int flatGlobal,
                                               const unordered_set<std::string>& levels); 

            //! Bias usually set to be zero
            vector< map<int, int> > etaBias_;

            //! Flag indicating MLWENOPrepare is used.
            bool prepare_ = false;
    };
}

#endif
