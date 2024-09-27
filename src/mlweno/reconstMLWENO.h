#ifndef RECONSTMLWENO_H_
#define RECONSTMLWENO_H_

#include "polynomial.h"
#include <map>

/**!
 * A new mlweno reconstuction header file.
 * An alternative to reconstruction.h file.
 * DONOT compile it with reconstruction.h. with cause fatal error.
 */

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

            double CalculateSmoothnessIndic(const MeshInfo& mi, indice owner, double** lu);
            double CalculateSmoothnessIndic(const MeshInfo& mi, indice owner, double** lu, const std::string& name);

            /**
             * Directly extract pre-calculateed smoothness indicator.
             * Should always be the one to call when smoothness indicator is needed.
             */
            double GetSmoothnessIndic(const MeshInfo& mi, indice owner); 

            // Extract pre-calculated smoothness indicator for system case
            double GetSmoothnessIndic(const MeshInfo& mi, indice owner, const std::string& name);

            unordered_map<int, double> GetSmoothnessIndicDeriv(const MeshInfo& mi, const indice& owner);

            //! Update smoothness indicator for entire reconstruction level
            void UpdateSmoothnessIndic(const MeshInfo& mi);

            //! Update smoothness indicator when transporting system of scalars
            void UpdateSmoothnessIndic(const MeshInfo& mi, double** lu, const std::string& name);

            void UpdateDerivSmoothnessIndic(const MeshInfo& mi);
            /**
             * Evaluate polynomial
             */
            double Evaluate(const MeshInfo& mi, const indice& owner, const vertex& point);

            double Evaluate(const MeshInfo& mi, const indice& owner, const vertex& point, const std::string& name);

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
            map<std::string,map<int, double>> smoothnessIndicVec_;          //!< Smoothenss indicators for system transport

            unordered_map<int, derivative> smoothnessIndicDeriv_;           //!< Derivative of smoothness indicators
            map<std::string, unordered_map<int, derivative>> smoothnessIndicDerivVec_;
                                                                            //!< Derivatives of smooth indic for system transport
            map<int, tensorProductPoly::stencilPolynomial*> singleLevel_;   //!< Single level polynomials

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

            void UpdateSmoothnessIndic(const MeshInfo& mi, double** lu, 
                                       const std::string& name);

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

            //! A destructor
            /**
             * Construct multi-level weno reconstruction by specifying each single level
             */
            ~multiLevelReconstruction() {};

            /**!
             * Select levels will be used.
             */
            void SelectWenoReconstLevel(const unordered_set<std::string>& keys,
                                        const MLWENOPrepare& mlpPtr);

            /**!
             * Assign reconstruction method to different levels.
             */
            void ModifyReconstMethod(const std::string& key,
                                     const vector<indice>& newReconstMethod);

            /**!
             * Assign linear weights to different levels.
             * This linear weights only used in new multilevel algorithm.
             */
            void SetUpLinearWgts(const std::string& key,
                                 const vector<double>& linwgts);

            /**!
             * Update Non linear weights.
             */
            void UpdateNonLinearWgts(const MeshInfo& mi,
                                     const std::string& weightType,
                                     bool (*assignML)(const indice& globalCell,
                                                      const MeshInfo& mi));

            /**!
             * Update Non linear weights.
             * No differentiating two-stage and one-stage
             */
            void UpdateNonLinearWgts(const MeshInfo& mi,
                                     bool (*assginML)(const indice& globalCell,
                                                      const MeshInfo& mi));

            /**!
             * Update Non linear weights.
             * Used for system transport.
             */

            void UpdateNonLinearWgts(const MeshInfo& mi,
                                     bool (*assginML)(const indice& globalCell,
                                                      const MeshInfo& mi),
                                     const std::string& name);

            /**!
             * Reconstruct point value with pre defined multi level WENO 
             * reconstruction scheme.
             */
            double EvaluateMLWENO (const MeshInfo& mi, 
                                   const vertex& point, 
                                   const indice& globalCell) const ;

            double EvaluateMLWENO (const MeshInfo& mi, 
                                   const vertex& point,
                                   const indice& globalCell,
                                   const std::string& name) const;

            /**!
             * Print information of non linear weights
             */
            void PrintNonLinearWgts(const MeshInfo& mi);

       private:

            double eps0_ = 1e-4;

            /**!
             * Storing all non linear weights mapping to each cells.
             * An efficient way to reuse calculated non linear weights.
             */
            unordered_map<int, unordered_map<std::string, unordered_map<int, double>>> nonLinearWgts_;

            unordered_map<std::string, unordered_map<int, unordered_map<std::string, unordered_map<int, double>>>> nonLinearWgtsVec_;

            /**!
             * Weno Levels and corresponding methods
             */
            unordered_map<std::string, singleLevelReconstruction *> Levels_;

            unordered_map<std::string, vector<indice>> Methods_;

            unordered_map<std::string, vector<double>> LinWgts_;

            /**!
             * Update non linear weight for one target cell.
             */
            void UpdateNonLinearWgtsCell_(const MeshInfo& mi,
                                          const int& globalCell,
                                          const std::string& weightType);

            /**!
             * Update nonlinear weights for one target cell.
             * Follow MLWENO paper.
             * Most recent definition of nonlinear weights.
             */
            void UpdateNonLinearWgtsCell_(const MeshInfo& mi, 
                                          const int& globalCell);

            /**!
             * Update nonlinear weights for the system transport
             */
            void UpdateNonLinearWgtsCell_(const MeshInfo& mi, 
                                          const std::string& name,
                                          const int& globalCell);
   };

}

#endif
