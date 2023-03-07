#ifndef SLRECONSTRUCTION_H_
#define SLRECONSTRUCTION_H_

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
            ~singleLevelReconstruction() {interior_.clear(); singleLevel_.clear();smoothnessIndic_.clear();};

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
            map<int, stencilPolynomial*> singleLevel_; //!< Single level polynomials
            map<int, double> smoothnessIndic_;         //!< Smoothness Indicators

            void IdentifyInteriorCell_(const MeshInfo& mi);

            const vertex ComputeStencilCenter_(const MeshInfo& mi, int flat);

            void ComputeStencilPolyn_(const MeshInfo& mi);

            void UpdateSmoothnessIndic_(const MeshInfo& mi);

    };

}

#endif
