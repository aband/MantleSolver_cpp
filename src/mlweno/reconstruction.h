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

namespace MLWENO {

    class reconstruction {

        public:
            reconstruction(){};
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
}
#endif
