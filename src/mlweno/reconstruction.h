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
#include <iomanip>

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

            void PrintStencils() const;

            void Clear();

        private:
            vector<int*> stencilSize_;
            
            void UpdateStencilSizeMax_(int* newStencilSize);

            int stencilSizeMax_[2] {0,0};
            vector<indice> shift_;

            int stencilNum_ = 0;

            vector<stencil <indice>> stencilIndice_;
            vector<stencilPolynomial*> stencilPolyn_;
    };

}
#endif
