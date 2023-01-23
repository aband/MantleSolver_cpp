#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

/*
 *A ML weno reconstruction contains several information:
 *1. Center point
 *2. Stencils
 *3. max polynomial order
 *
 *What happens in the class of Reconstruction:
 *1. Compute non linear weights
 *2. Compute reconstruction values
 *3. Compute derivatives of reconstruction values
 */

#include "polynomial.h"
#include <map>

namespace MLWENO {

    class reconstruction {

        public:
            reconstruction(){};
            // Create stencil polynomials from stencil of indices
            reconstruction(vector<stencil <indice>> stencilIndice) {};
            reconstruction(int* stencilSize, indice shift) {AddStencil(stencilSize, shift);};
            reconstruction(vector<int*> stencilSizes, vector<indice> shifts) {AddStencil(stencilSizes,shifts);};

            ~reconstruction(){};

            void AddStencil(int* stencilSize, indice shift);
            void AddStencil(vector<int*> stencilSizes, vector<indice> shifts);

            void CreateStencilPolynomials(const indice& start,              const vertex& center,
                                          const vector<indice>& targetCell, const MeshInfo& mi);

            void AddStencilPolynomials(const indice& start,              const vertex& center,
                                       const vector<indice>& targetCell, const MeshInfo& mi);

            void Update(const MeshInfo& mi) {ComputeNonLinWgts_(mi);};

            double Eval(double x, double y) const;
            double Eval(const vertex& P) const {return Eval(P[0],P[1]);};
            double operator() (const double x, const double y) const {return Eval(x,y);};
            double operator() (const vertex& P) const {return Eval(P);};

            // Print out calculated private variables
            void PrintSmoothnessIndic() const;
            void PrintStencils() const;
            void PrintNonLinWgts() const;

            // Clear class member vectors and reset to default values
            void Clear();

        private:

            double eps0_ = 0.01;
 
            double scale_ = -1;

            vector<int> etaBias_;

            vector<int*> stencilSize_;
            
            void UpdateStencilSizeMax_(int* newStencilSize);

            int stencilSizeMax_[2] {0,0};
            vector<indice> shift_;

            int stencilNum_ = 0;

            vector<stencil <indice>> stencilIndice_;
            vector<stencilPolynomial*> stencilPolyn_;

            void ComputeNonLinWgts_(const MeshInfo& mi);
            void ComputeSmoothnessIndicatorPolyn_(const MeshInfo& mi);

            vector<double> linWgts_;
            vector<double> nonLinWgts_;

            vector<double> smoothnessIndicPolyn_;
    };

    // A reconstruction methodology based on the idea of multi level ideas

    class singleLevelReconstruction {
        public:
            singleLevelReconstruction() {};
            singleLevelReconstruction(int stencilSizeX, int stencilSizeY); 
            
            ~singleLevelReconstruction() {interior_.clear(); singleLevel_.clear();smoothnessIndic_.clear();};

            void CreateStencilPolynomials(const MeshInfo& mi);

            void CheckStencils() const {cout<< "Constructed "<< interior_.size() << " stencils with the size of " << stencilSizeX_ << " " << stencilSizeY_ << endl;};
            void CheckStencilPolynomials(const MeshInfo& mi, indice start);

            const double CalculateSmoothnessIndic(const MeshInfo& mi, indice owner);

            const double GetScale(const MeshInfo& mi, indice owner) {return singleLevel_[FlatIndic(mi,owner)]->GetScale();}; 

            int CheckExist(const MeshInfo& mi, indice owner) const {return interior_.count(FlatIndic(mi,owner));};

        private:

            vertex ComputeStencilCenter_(const MeshInfo& mi, int flat);

            int stencilSizeX_ = -1;
            int stencilSizeY_ = -1;

            void IdentifyInteriorCell_(const MeshInfo& mi);
            unordered_set<int> interior_;

            void ComputeStencilPolyn_(const MeshInfo& mi);
            map<int, stencilPolynomial*> singleLevel_;

            void UpdateSmoothnessIndic_(const MeshInfo& mi);
            stencil <indice> stencilIndice_;

            map<int, double> smoothnessIndic_;

    };

    class multiLevelReconstruction {
        public:
            multiLevelReconstruction() {};
            multiLevelReconstruction(const MeshInfo& mi, int stencilSizeX, int stencilSizeY){AddLevel(mi,stencilSizeX,stencilSizeY);};
            multiLevelReconstruction(const MeshInfo& mi, int* stencilSize) {AddLevel(mi,stencilSize);};
            multiLevelReconstruction(const MeshInfo& mi, vector<int*> stencilSizes) 
            {AddLevel(mi,stencilSizes);};

            ~multiLevelReconstruction() {Clear();};

            void AddLevel(const MeshInfo& mi, int stencilSizeX, int stencilSizeY);
            void AddLevel(const MeshInfo& mi, int* stencilSize) {AddLevel(mi,stencilSize[0],stencilSize[1]);};
            void AddLevel(const MeshInfo& mi, vector<int*> stencilSizes)
            {for (int i=0; i<stencilSizes.size(); i++){
                 AddLevel(mi,stencilSizes[i]);}};

            void AddReconstMethod(vector<indice> brm) {baseReconstMethod_.push_back(brm); AddWgts_();}; 
            void AddReconstMethod(vector<vector<indice>> brms) 
                                 {for (int i=0; i<brms.size(); i++){
                                      AddReconstMethod(brms[i]);     
                                  }}; 

            void UpdateNonlinearWgts(const MeshInfo& mi, indice start);

            void GetInfo();

            void Clear();
        private:

            const double eps0_ = 0.001;

            vector< singleLevelReconstruction *> allLevels_;

            vector<vector<indice>> baseReconstMethod_; 

            void ResetWgts_();
            void AddWgts_();

            vector< vector<double> > linearWgts_;
            map<int, vector< vector<double> >> nonLinearWgts_;

            vector< vector<int> > etaBias_;
    };
}
#endif
