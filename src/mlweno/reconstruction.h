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
            
            ~singleLevelReconstruction() {interior_.clear(); singleLevel_.clear();};

            void CreateStencilPolynomials(const MeshInfo& mi);

            void CheckStencils() const {cout<< "Constructed "<< interior_.size() << " stencils with the size of " << stencilSizeX_ << " " << stencilSizeY_ << endl;};
            void CheckStencilPolynomials(const MeshInfo& mi, indice start);

        private:

            vertex ComputeStencilCenter_(const MeshInfo& mi, int flat);

            // Flatten indice into 1D array
            int FlatIndic_(const MeshInfo& mi, int i, int j) const 
                          {return j*mi.MPIlocalCellSize[0]+i;};

            int FlatIndic_(const int M, int i, int j) const {return j*M+i;};
            int FlatIndic_(const MeshInfo& mi, const indice& p) const {return FlatIndic_(mi,p[0],p[1]);}
            int FlatIndic_(const int M, const indice& p) const {return FlatIndic_(M,p[0],p[1]);};

            // Reverse process of flatten indices
            indice Bend_(const MeshInfo& mi, int flat) const 
                        {return {flat%mi.MPIlocalCellSize[0], flat/mi.MPIlocalCellSize[0]};};

            indice Bend_(const int M, int flat) const {return {flat%M, flat/M};}

            int stencilSizeX_ = -1;
            int stencilSizeY_ = -1;

            void IdentifyInteriorCell_(const MeshInfo& mi);
            unordered_set<int> interior_;

            void ComputeStencilPolyn_(const MeshInfo& mi);
            map<int, stencilPolynomial*> singleLevel_;

            stencil <indice> stencilIndice_;

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

            
            


            void UpdateNonlinearWgts();

            void GetInfo();

            void Clear();
        private:

            vector< singleLevelReconstruction *> allLevels_;

            vector<vector<indice>> reconstMethod_; 

            vector< vector<double> > linearWgts_;
            vector< vector<double> > nonLinearWgts_;
    };
}
#endif
