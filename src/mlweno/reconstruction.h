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

#include <vector>
#include <pair>
#include <set>
#include <unordered_set>
#include <array>
#include <valarray>
#include <algorithm>
#include <numeric>
#include <memory>

#include "lapacke.h"
#include "integral.h"
#include <assert.h>


namespace MLWENO {

    class reconstruction {

        public:
            reconstruction(){};
            reconstruction(int[2] stencilSize, int[2] shift) {AddStencil(stencilSize, shift);};
            ~reconstruction(){};

            void AddStencil(int[2] stencilSize, int[2] shift);
            void AddStencil(vector<int[2]> stencilSizes, vector<int[2]> shifts);

        private:
            vector<int[2]> stencilSize_;
            
            void UpdateStencilSizeMax_(int[2] newStencilSize);

            int stencilSizeMax_[2] {0,0};
            vector<int[2]> shift_;

            int stencilNum_ = 0;

            vector<stencil <indice>> * stenilIndice_ = nullptr;
            vector<stencilPolynomial> stncilPolyn_;
    }

}
#endif
